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

#include "LinearHeatIntegrator.h"
#include "amsiFEA.h"
#include "gmi.h"

namespace mumfim
{
  LinearHeatConductionStep::LinearHeatConductionStep(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm comm_)
      : amsi::FEAStep(mesh, analysis_case, {}, {}, "macro", comm_)
      , constitutives()
      , iteration(0)
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh, "temperature", apf::SCALAR, 1);
    apf::zeroField(apf_primary_field);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);

    kappa = apf::createIPField(apf_mesh, "kappa", apf::MATRIX, 1);
    apf::zeroField(kappa);
    // used to place kappa into IPfield by region
    std::map<int, apf::Matrix3x3> mappa;

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
      const auto * kappa_mt =
          //mt::GetCategoryModelTraitByType<mt::MatrixMT>(continuum_model,
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "kappa");
      if (kappa_mt == nullptr)
      {
        std::cerr << " \"kappa\" (thermal conductivity) is required for "
                     "the continuum model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }

      auto k = (*kappa_mt)();
      // This needs to be unique_ptr-ized, as it currently isn't deleted anywhere.
      apf::Matrix3x3 * kappa_r = new apf::Matrix3x3(
        k  , 0.0, 0.0,
        0.0, k  , 0.0,
        0.0, 0.0, k  );
      // Save for the kappa IPField assignment
      mappa[tag] = *kappa_r;

      std::cout << "continuum model type: " << continuum_model->GetType()
                << "\n";
      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<LinearHeatIntegrator>(
              apf_primary_field, apf_primary_numbering, kappa_r);
    }
    gmi_end(gmodel, it);

    apf::MeshEntity *ent;
    auto *mesh_it = mesh->begin(3);
    while((ent = mesh->iterate(mesh_it)))
    {
      int tag = mesh->getModelTag(mesh->toModel(ent));
      apf::setMatrix(kappa, ent, 0, mappa.at(tag));
    }

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
    AssembleIntegratorIntoLAS(las);
  }


  /* Linear problem, just solve it as the initial guess*/
  void LinearHeatConductionStep::computeInitGuess(amsi::LAS * las)
  {
    setSimulationTime(T);
    LinearSolver(this, las);
    apf::writeASCIIVtkFiles("foo", this->apf_mesh);
    las->iter();
  }

  amsi::ElementalSystem * LinearHeatConductionStep::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
    return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }
}