#include "LinearTissueStep.h"

#include <amsiLinearElasticConstitutive.h>
#include <apfFunctions.h>
#include <gmi.h>

#include "amsiFEA.h"

namespace mumfim
{
  LinearTissueStep::LinearTissueStep(apf::Mesh * mesh,
                            const mt::CategoryNode & analysis_case,
                            MPI_Comm cm)
      : amsi::FEAStep(mesh,
                   analysis_case,
                   {},
                   {},
                   "macro",
                   cm)
      , constitutives()
  {
    apf_primary_field = apf::createLagrangeField(
        apf_mesh, "linear_displacement", apf::VECTOR, 1);
    apf::zeroField(apf_primary_field);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    // FIXME Dirichlet BC entry should take a mt_type "displacement" to make
    // this explicit
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
    auto * gmodel = mesh->getModel();
    if (gmodel == nullptr)
    {
      std::cerr << "Mesh must have attached model!\n";
      exit(1);
    }
    static constexpr int dimension = 3;
    auto * it = gmi_begin(gmodel, dimension);
    struct gmi_ent * gent;
    while ((gent = gmi_next(gmodel, it)))
    {
      int tag = gmi_tag(gmodel, gent);
      const auto * region_traits =
          problem_definition.associated->Find({dimension, tag});
      if (region_traits == nullptr)
      {
        std::cerr << "region " << tag << "must have material model data\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * material_model =
          mt::GetCategoryByType(region_traits, "material model");
      if (material_model == nullptr)
      {
        std::cerr << "material model must exist on region" << tag << "\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * continuum_model =
          mt::GetPrimaryCategoryByType(material_model, "continuum model");
      if (material_model == nullptr)
      {
        std::cerr << "continuum model must exist on region" << tag << "\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * youngs_modulus =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "youngs modulus");
      const auto * poisson_ratio =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "poisson ratio");
      if (youngs_modulus == nullptr || poisson_ratio == nullptr)
      {
        std::cerr << "youngs modulus or poisson ratio missing on region " << tag
                  << "\n";
        exit(1);
      }
      constitutives[tag] = std::make_unique<amsi::LinearElasticIntegrator>(
          apf_primary_field, apf_primary_numbering, 1, (*youngs_modulus)(), (*poisson_ratio)());
    }
    gmi_end(gmodel, it);
  }
  void LinearTissueStep::UpdateDOFs(const double * sl)
  {
    amsi::WriteOp wrop;
    amsi::FreeApplyOp frop(apf_primary_numbering, &wrop);
    amsi::ApplyVector(apf_primary_numbering, apf_primary_field, sl,
                      first_local_dof, &frop)
        .run();
    apf::synchronize(apf_primary_field);
  }
  void LinearTissueStep::Assemble(amsi::LAS * las)
  {
    // For the LinearTissueStep, the coordinate field is assumed to be the meshes
    // coordinate field. This is unlike what's used in NonlinearTissueStep and
    // MultiscaleTissueStep
    // Note that AssembleIntegratorIntoLAS zeros out the LAS matrix 
    // and vector, so add in the NeumannBCs after, not before.
    AssembleIntegratorIntoLAS(las);
    ApplyBC_Neumann(las);
  }

  amsi::ElementalSystem * LinearTissueStep::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
    return constitutives[apf_mesh->getModelTag(apf_mesh->toModel(mesh_entity))]
        .get();
  }
}  // namespace mumfim
