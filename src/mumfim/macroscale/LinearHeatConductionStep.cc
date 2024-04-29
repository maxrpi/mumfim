#include "LinearHeatConductionStep.h"

#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfLabelRegions.h>

#include <array>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gmi.h"
#include "LinearHeatIntegrator.h"

namespace mumfim
{
  LinearHeatConductionStep::LinearHeatConductionStep(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm comm_)
      : AnalysisStep(mesh, analysis_case, {}, {}, "macro", comm_)
      , constitutives()
      , iteration(0)
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh, "temperature", apf::SCALAR, 1);
    apf::zeroField(apf_primary_field);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);

    kappa = apf::createIPField(apf_mesh, "kappa", apf::MATRIX, 1);

    // Future expansion
    //flux = apf::createIPField(apf_mesh, "fluxes", apf::VECTOR, 1);
    //sources = apf::createIPField(apf_mesh, "sources", apf::SCALAR, 1);
    //apf::zeroField(sources);
    //apf::zeroField(flux);

    amsi::applyUniqueRegionTags(apf_mesh);
    
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
      const auto * tmp =
          //mt::GetCategoryModelTraitByType<mt::MatrixMT>(continuum_model,
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "kappa");
      auto k = (*tmp)();
      // This needs to be unique_ptr-ized, as it currently isn't deleted anywhere.
      apf::Matrix3x3 * kappa_r = new apf::Matrix3x3(
        k  , 0.0, 0.0,
        0.0, k  , 0.0,
        0.0, 0.0, k  );

      /*
      apf::Matrix3x3 * kappa_r = new apf::Matrix3x3(
        (*tmp)(0,0), (*tmp)(1,0), (*tmp)(2,0),
        (*tmp)(0,1), (*tmp)(1,1), (*tmp)(2,1),
        (*tmp)(0,2), (*tmp)(1,2), (*tmp)(2,2)
      );
      */

      std::cout << "continuum model type: " << continuum_model->GetType()
                << "\n";
      if (kappa == nullptr)
      {
        std::cerr << " \"kappa\" (thermal conductivity) is required for "
                     "the continuum model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<LinearHeatIntegrator>(apf_primary_field, kappa_r);
    }
    gmi_end(gmodel, it);
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{ .categories = {"heat_flux"},
                              .mt_name = "x",
                              .mt_type = amsi::NeumannBCType::pressure});
                              //.mt_type = amsi::NeumannBCType::normal_flux});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{ .categories = {"heat_flux"},
                              .mt_name = "y",
                              .mt_type = amsi::NeumannBCType::pressure});
                              //.mt_type = amsi::NeumannBCType::normal_flux});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{ .categories = {"heat_flux"},
                              .mt_name = "z",
                              .mt_type = amsi::NeumannBCType::pressure});
                              //.mt_type = amsi::NeumannBCType::normal_flux});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{ .categories = {"temperature"},
                                .mt_name = "temperature"});
  }

  LinearHeatConductionStep::~LinearHeatConductionStep()
  {
    apf::destroyField(apf_primary_field);
    apf::destroyField(flux);
    apf::destroyField(kappa);
    apf::destroyField(sources);
  }
  void LinearHeatConductionStep::Assemble(amsi::LAS * las)
  {
    ApplyBC_Neumann(las);
    AssembleIntegratorIntoLAS(
        las,
        [this](apf::MeshEntity * me, int) {
          auto tmp = constitutives[apf_mesh->toModel(me)].get();
          return tmp;
        });
  }


  /* Linear problem, just solve it as the initial guess*/
  void LinearHeatConductionStep::computeInitGuess(amsi::LAS * las)
  {
    setSimulationTime(T);
    LinearSolver(this, las);
    apf::writeASCIIVtkFiles("foo", this->apf_mesh);
    las->iter();
  }
  /*
  void LinearHeatConductionStep::step()
  {
    iteration = 0;
  }
  */
  /*
  void LinearHeatConductionStep::iter()
  {
    iteration++;
  }
  */

  //void LinearHeatConductionStep::UpdateDOFs(const double * sol) { ; }

  /**
   * get the Neumann boundary condition value on the specified entity
   * @param ent model entity to query
   * @param frc will be filled with the value of the NeumannBCs
   * This really needs to be reimplemented to be heat flux specific, but I
   * wasn't willing to just throw away the functionality without replacing it.
   */
  /*
  void LinearHeatConductionStep::getFluxOn(apf::ModelEntity * ent, double * frc)
  {
    auto model_dimension = apf_mesh->getModelType(ent);
    auto model_tag = apf_mesh->getModelTag(ent);
    auto * associated_traits =
        problem_definition.associated->Find({model_dimension, model_tag});
    if (associated_traits == nullptr)
    {
      std::cerr << "Attempting to get Neumann BC on [" << model_dimension << ","
                << model_tag << "].";
      std::cerr << " No boundary conditions associated with that geometry.\n";
      exit(1);
    }
    for (const auto & neumann_bc : neumann_bcs)
    {
      const mt::AssociatedCategoryNode * category_node = associated_traits;
      for (const auto & category : neumann_bc.categories)
      {
        category_node = category_node->FindCategoryByType(category);
        if (category_node == nullptr)
        {
          break;
        }
      }
      if (category_node == nullptr)
      {
        std::cerr << "Invalid Neumann BC provided\n";
        exit(1);
      }
      const auto * flux_trait =
          mt::GetCategoryModelTraitByType(category_node, neumann_bc.mt_name);
      const auto * const_flux_trait = mt::MTCast<mt::VectorMT>(flux_trait);
      const auto * time_function_flux =
          mt::MTCast<mt::VectorFunctionMT<1>>(flux_trait);
      if (const_flux_trait != nullptr)
      {
        for (int i = 0; i < 3; ++i)
        {
          frc[i] = (*const_flux_trait)(i);
        }
      }
      else if (time_function_flux != nullptr)
      {
        for (int i = 0; i < 3; ++i)
        {
          frc[i] = (*time_function_flux)(i, T);
        }
      }
      else
      {
        std::cerr << "flux must be either a scalar value, or function of time "
                     "only.\n";
        exit(1);
      }
    }
  }
*/
/*
  void LinearHeatConductionStep::recoverSecondaryVariables(int)
  {
    ;
  }
  */

}  // namespace mumfim