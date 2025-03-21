#include "LinearHeatConductionStep.h"

#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfFieldData.h>
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

    bool dumpqpmap = false, readqmap = false;
    if(dumpqpmap)
    {
      FILE *fp;
      fp = fopen("quadpointmap.txt", "w");

      apf::MeshTag* mesh_tag = apf_mesh->createIntTag("index", 1);
      int index[1]= {0};
      apf::MeshEntity *ent;
      auto *mesh_it = apf_mesh->begin(3);
      while((ent = apf_mesh->iterate(mesh_it)))
      {
        apf_mesh->setIntTag(ent, mesh_tag, index);
        int model_tag = apf_mesh->getModelTag(mesh->toModel(ent));
        apf::MeshElement* e = apf::createMeshElement(apf_mesh,ent);
        apf::Vector3 point;
        apf::getIntPoint(e,0,0,point);
        apf::Vector3 xyz;
        apf::mapLocalToGlobal(e, point, xyz);
        apf::destroyMeshElement(e);
        fprintf(fp, "%g, %g, %g, %d, %d\n", xyz[0], xyz[1], xyz[2], model_tag, *index);
        index[0]++;
      }
      fclose(fp);
    } 
    if(readqmap) 
    {
      FILE *fp;
      fp = fopen("kappavalues.txt", "r");
      int mesh_index;
      std::vector<apf::Matrix3x3*> kappa_me; kappa_me.reserve(apf_mesh->count(3));
      double kxx, kyy, kzz, kxy, kyz, kzx;
      while(fscanf(fp, "%d, %lf, %lf, %lf, %lf, %lf, %lf\n",
            &mesh_index, &kxx, &kyy, &kzz, &kxy, &kyz, &kzx) == 7)
      {
        apf::Matrix3x3 *kmatrix = new apf::Matrix3x3(
          kxx, kxy, kzx,
          kxy, kyy, kyz,
          kzx, kyz, kzz
        );
        kappa_me[mesh_index] = kmatrix; 
      }
      fclose(fp);
      mesh_index = 0;
      apf::MeshEntity *ent;
      auto *mesh_it = apf_mesh->begin(3);
      while((ent = apf_mesh->iterate(mesh_it)))
      {
        apf::setMatrix(kappa, ent, 0, *kappa_me[mesh_index]);
        mesh_index++;
      }
      kappa_me.clear();
    } 

    // used to place kappa into IPfield by region
    std::map<int, apf::Matrix3x3> mappa;

    amsi::applyUniqueRegionTags(apf_mesh);
    
    static constexpr int dimension = 3;
    auto * gmodel = apf_mesh->getModel();
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
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "kappa");
      const auto * kappa_tensor_mt =
          mt::GetCategoryModelTraitByType<mt::VectorMT>(continuum_model,
                                                        "kappaTensor");

      if(kappa_mt != nullptr && kappa_tensor_mt != nullptr){
        std::cerr << "Specifying both \"kappa\" and \"kappaTensor\" for the "
                      "same material is not allowed.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      if (!readqmap && kappa_mt == nullptr && kappa_tensor_mt == nullptr)
      {
        std::cerr << " \"kappa\" or \"kappaTensor\" (thermal conductivity) is "
                     "required for the continuum model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }

      // This needs to be unique_ptr-ized, as it currently isn't deleted anywhere.
      apf::Matrix3x3 * kappa_r;

      if(kappa_mt != nullptr){
        auto k = (*kappa_mt)();
        kappa_r = new apf::Matrix3x3(
          k  , 0.0, 0.0,
          0.0, k  , 0.0,
          0.0, 0.0, k  );
      } else {
        kappa_r = new apf::Matrix3x3(
          (*kappa_tensor_mt)(0), (*kappa_tensor_mt)(1), (*kappa_tensor_mt)(2), 
          (*kappa_tensor_mt)(3), (*kappa_tensor_mt)(4), (*kappa_tensor_mt)(5), 
          (*kappa_tensor_mt)(6), (*kappa_tensor_mt)(7), (*kappa_tensor_mt)(8)
        );
      }
      // Save for the kappa IPField assignment
      mappa[tag] = *kappa_r;

      if(readqmap)
      {
      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<LinearHeatIntegrator>(
              apf_primary_field, apf_primary_numbering, kappa);
      } else 
      {
      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<LinearHeatIntegrator>(
              apf_primary_field, apf_primary_numbering, kappa_r);
      }
    }
    gmi_end(gmodel, it);

    /*
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(3);
    while((ent = apf_mesh->iterate(mesh_it)))
    {
      int model_tag = apf_mesh->getModelTag(mesh->toModel(ent));
      apf::setMatrix(kappa, ent, 0, mappa.at(model_tag));
    }
    */

    
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{ .categories = {"heat_flux"},
                              .mt_name = "flux",
                              .mt_type = amsi::NeumannBCType::pressure});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{ .categories = {"convection"},
                              .mt_name = "parameters",
                              .mt_type = amsi::NeumannBCType::robin});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{ .categories = {"temperature"},
                                .mt_name = "T"});
  }

  LinearHeatConductionStep::~LinearHeatConductionStep()
  {
    apf::destroyField(apf_primary_field);
    apf::destroyField(kappa);
  }
  void LinearHeatConductionStep::Assemble(amsi::LAS * las)
  {
    // Note that AssembleIntegratorIntoLAS zeros out the LAS matrix 
    // and vector, so add in the NeumannBCs after, not before.
    AssembleIntegratorIntoLAS(las);
    ApplyBC_Neumann(las);
  }

  void LinearHeatConductionStep::UpdateDOFs(const double * solution)
  {
    int num_components = apf::countComponents(apf_primary_field);
    assert(num_components == 1 && "Nonscalar field for temperature");

    apf::MeshEntity * mesh_ent = NULL;
    for (int ii = 0; ii < analysis_dim; ii++)
    {
      for (apf::MeshIterator * it = apf_mesh->begin(ii);
           (mesh_ent = apf_mesh->iterate(it));) {

        if (apf_mesh->isOwned(mesh_ent)) {
          apf::FieldShape * shape_function = apf::getShape(apf_primary_field);

          for (int jj = 0;
               jj < shape_function->countNodesOn(apf_mesh->getType(mesh_ent));
               jj++) {

            if (!apf::isFixed(apf_primary_numbering, mesh_ent, jj, 0)) {
              int global_number = getNumber(apf_primary_numbering, mesh_ent, jj, 0);
              auto temperature = solution[global_number - first_local_dof];
              apf::setScalar(apf_primary_field, mesh_ent, jj, temperature);
            }
          }
        }
      }
    }
    apf::synchronize(apf_primary_field);
  }


  /* Linear problem, just solve it as the initial guess*/
  void LinearHeatConductionStep::computeInitGuess(amsi::LAS * las)
  {
    setSimulationTime(T);
    LinearSolver(this, las);
    las->iter();
  }

  amsi::ElementalSystem * LinearHeatConductionStep::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
    return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }
}