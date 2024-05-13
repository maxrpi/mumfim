#include "NonlinearTissueStep.h"

#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfLabelRegions.h>

#include <array>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include "NeoHookeanIntegrator.h"
#include "amsiFEA.h"
#include "gmi.h"

namespace mumfim
{
  NonlinearTissueStep::NonlinearTissueStep(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm cm)
      : amsi::FEAStep(mesh, analysis_case, {}, {}, "macro", cm)
      , constitutives()
      , dv_prev(0.0)
      , load_step(0)
      , iteration(0)
      , _analysis_case(analysis_case)
  {
    apf_primary_field =
        apf::createLagrangeField(apf_mesh, "displacement", apf::VECTOR, 1);
    apf::zeroField(apf_primary_field);
    delta_u = apf::createLagrangeField(apf_mesh, "displacement_delta",
                                       apf::VECTOR, 1);
    apf::zeroField(delta_u);
    apf_primary_delta_field = delta_u;

    accepted_displacements = apf::createLagrangeField(
        apf_mesh, "accepted_displacement", apf::VECTOR, 1);
    apf::zeroField(accepted_displacements);
    // create a current coordinate field from the CurrentCoordFunc (x = X+u)
    xpyfnc =
        new amsi::XpYFunc(apf_mesh->getCoordinateField(), apf_primary_field);
    // user field is not zeroed because it is the sum of two other fields
    current_coords =
        apf::createUserField(apf_mesh, "current_coordinates", apf::VECTOR,
                             apf::getLagrange(1), xpyfnc);
    // take the current coordinate field, and subtract increment in displacement
    // to get the coords at the previous increment (this is sorta hacky but it
    // should work)
    prv_crd_fnc = new amsi::XpYFunc(current_coords, delta_u, 1, -1);
    prev_coords =
        apf::createUserField(apf_mesh, "previous_coordinates", apf::VECTOR,
                             apf::getLagrange(1), prv_crd_fnc);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    strs = apf::createIPField(apf_mesh, "stress", apf::MATRIX, 1);
    apf::zeroField(strs);
    strn = apf::createIPField(apf_mesh, "strain", apf::MATRIX, 1);
    apf::zeroField(strn);
    dfm_grd = apf::createIPField(apf_mesh, "F", apf::MATRIX, 1);
    apf::zeroField(dfm_grd);
    stf_vrtn =
        apf::createIPField(apf_mesh, "stiffness_variation", apf::SCALAR, 1);
    apf::zeroField(stf_vrtn);
    axl_yngs_mod =
        apf::createIPField(apf_mesh, "axial_youngs_modulus", apf::SCALAR, 1);
    apf::zeroField(axl_yngs_mod);
    amsi::applyUniqueRegionTags(apf_mesh);
    // zero stiffness_variation and axial_youngs_modulus fields
    const auto * stiffness_gradient =
        problem_definition.unassociated->FindCategoryNodeByType(
            "stiffness gradient");
    if (stiffness_gradient != nullptr)
    {
      for (const auto & grad_nd : stiffness_gradient->GetCategoryNodes())
      {
        stf_vrtn_cnst.push_back(buildStiffnessVariation(grad_nd, stf_vrtn));
      }
    }
    for (auto & cnst : stf_vrtn_cnst)
    {
      cnst->populate_stf_vrtn_fld();
    }
    const auto * incompressible =
        problem_definition.unassociated->FindCategoryNodeByType(
            "incompressible");
    if (incompressible != nullptr)
    {
      throw mumfim_error("incompressibility not currently supported");
    }
    static constexpr int dimension = 3;
    auto * gmodel = mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, dimension);
    while ((gent = gmi_next(gmodel, it)))
    {
      int tag = gmi_tag(gmodel, gent);
      const auto * model_traits =
          problem_definition.associated->Find({dimension, tag});
      const auto * material_model =
          mt::GetCategoryByType(model_traits, "material model");
      const auto * continuum_model =
          mt::GetPrimaryCategoryByType(material_model, "continuum model");
      const auto * youngs_modulus =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "youngs modulus");
      const auto * poisson_ratio =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "poisson ratio");
      // Even pytorch needs a young's modulus, poisson ratio b/c
      // LinearTissueStep is used as an initial guess. The initial guess should
      // be based off the same continuum material model as the nonlinear tissue
      // TODO: LinearTissueStep uses same continuum material model as
      // NonlinearTissueStep
      std::cout << "continuum model type: " << continuum_model->GetType()
                << "\n";
      if ((youngs_modulus == nullptr || poisson_ratio == nullptr))
      {
        std::cerr << " \"poisson ratio\" and \"youngs modulus\" are required "
                     "types for the continuum model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      if (continuum_model->GetType() == "isotropic_neohookean")
      {
        // should check to make sure the continuum model is iso lin ela for init
        // solve?
        constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
            std::make_unique<NeoHookeanIntegrator>(
                apf_primary_field, dfm_grd, current_coords, strs, strn,
                (*youngs_modulus)(), (*poisson_ratio)(), 1);
      }
      else if (continuum_model->GetType() == "transverse_isotropic")
      {
        throw mumfim_error("Transverse isotropic material not currently supported");
      }
      else if (continuum_model->GetType() == "pytorch")
      {
#ifdef MUMFIM_ENABLE_TORCH
        const auto * model_file = mt::GetCategoryModelTraitByType<mt::StringMT>(
            continuum_model, "model file");
        constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
            std::make_unique<
                UpdatedLagrangianMaterialIntegrator<TorchMaterial>>(
                TorchMaterial{(*model_file)()}, strn, strs, apf_primary_field,
                dfm_grd, 1);
#else
        throw mumfim_error("MUMFIM was not compiled with PyTorch support");
#endif
      }
      else {
        throw mumfim_error("Invalid continuum model was supplied");
      }
    }
    gmi_end(gmodel, it);
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "x"});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "y"});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "z"});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{.categories = {"pressure"},
                             .mt_name = "magnitude",
                             .mt_type = amsi::NeumannBCType::pressure});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{.categories = {"traction"},
                             .mt_name = "direction",
                             .mt_type = amsi::NeumannBCType::traction});
  }

  NonlinearTissueStep::~NonlinearTissueStep()
  {
    delete xpyfnc;
    apf::destroyField(current_coords);
    delete prv_crd_fnc;
    apf::destroyField(prev_coords);
    // destroying this field currently throws an error
    // most likely because it is zeroed and core does not
    // properly deal with this...
    // apf::destroyField(delta_u);
    apf::destroyField(strs);
    apf::destroyField(strn);
    apf::destroyField(stf_vrtn);
    apf::destroyField(axl_yngs_mod);
  }

  void NonlinearTissueStep::computeInitGuess(amsi::LAS * las)
  {
    // LinearTissueStep lt(model, mesh, prob_def, solution_strategy, analysis_comm);
    LinearTissueStep lt(apf_mesh, _analysis_case, analysis_comm);
    lt.setSimulationTime(T);
    LinearSolver(&lt, las);
    las->iter();
    apf::copyData(delta_u, lt.getUField());
    apf::copyData(apf_primary_field, lt.getUField());
  }

  void NonlinearTissueStep::step()
  {
    iteration = 0;
    load_step++;
  }

  void NonlinearTissueStep::iter()
  {
    iteration++;
  }

  void NonlinearTissueStep::Assemble(amsi::LAS * las)
  {
    ApplyBC_Neumann(las);
    // custom iterator would be perfect for switching for multiscale version
    // FIXME !!! here createMeshELement is using the mesh (Original coords),
    AssembleIntegratorIntoLAS(las, current_coords);
  }

  void NonlinearTissueStep::UpdateDOFs(const double * sol)
  {
    // rewrite for SNES. Here the solution vector we have is the
    // full displacements
    // FIXME: Detangle this code once SNES is working
    amsi::WriteOp wr_op;
    auto num_dofs = countNodes(apf_primary_numbering) *
                    countComponents(apf_primary_delta_field);

    std::vector<double> old_solution(num_dofs);
    amsi::ToArray(apf_primary_numbering, accepted_displacements,
                  &old_solution[0], first_local_dof, &wr_op)
        .run();
    amsi::FreeApplyOp frwr_op(apf_primary_numbering, &wr_op);
    amsi::ApplyVector(apf_primary_numbering, apf_primary_field, sol,
                      first_local_dof, &frwr_op)
        .run();
    std::transform(
        old_solution.begin(), old_solution.end(), sol, old_solution.begin(),
        [](double old_sol, double new_sol) { return new_sol - old_sol; });

    amsi::ApplyVector(apf_primary_numbering, apf_primary_delta_field,
                      old_solution.data(), first_local_dof, &frwr_op)
        .run();

    apf::synchronize(apf_primary_field);
    apf::synchronize(delta_u);
  }

  amsi::ElementalSystem * NonlinearTissueStep::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
        return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }

}  // namespace mumfim
