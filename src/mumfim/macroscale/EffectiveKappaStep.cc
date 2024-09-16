#include "EffectiveKappaStep.h"

#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfLabelRegions.h>

#include <array>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include "EffectiveKappaIntegrator.h"
#include "amsiFEA.h"
#include "gmi.h"

namespace mumfim
{
  EffectiveKappaStep::EffectiveKappaStep(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm comm_)
      : amsi::FEAStep(mesh, analysis_case, {}, {}, "macro", comm_)
      , constitutives()
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh, "dummy", apf::SCALAR, 1);
    apf::zeroField(apf_primary_field);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);

    kappa = apf::createIPField(apf_mesh, "kappa", apf::MATRIX, 1);
    apf::zeroField(kappa);
    
    // Develop a vector of model faces on the boundary
    std::vector<struct gmi_ent *> boundary_model_faces;
    auto * gmodel = apf_mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, 2); // Loop over model faces
    while ((gent = gmi_next(gmodel, it)))
    {
      gmi_set adjacent_regions = gmi_adjacent(gmodel, gent, 3);
      if(adjacent_regions->n == 1){
        boundary_model_faces.push_back(gent);
      }
    }
    gmi_end(gmodel, it);

    // Populate vectors of interior and exterior nodes
    apf::MeshEntity *ent;
    auto *mesh_it = mesh->begin(0);
    while((ent = mesh->iterate(mesh_it)))
    {
      apf::ModelEntity * classifiedOn = apf_mesh->toModel(ent);
      found = false;
      for(auto bmf = boundary_model_faces.begin();
          bmf != boundary_model_faces.end(); 
          bmf++)
      {
        if(apf::isInClosureOf(classifiedOn, bmf))
        {
          exterior.pushback(ent);
          found = true;
          break;
        }
      }
      if(!found){
        interior.push_back(ent);
      }
    }

    // Allocate submatrices of partitioned stiffness matrix
    int n_int = interior.size;
    int n_ext = exterior.size;
    assert(n_int + n_ext == apf_mesh->count(0));
    //These need to be a better matrix data structure. This are biggish and sparse.
    kii = apf::DynamicMatrix(n_int, n_int);
    kei = apf::DynamicMatrix(n_ext, n_int); 
    kie = apf::DynamicMatrix(n_int, n_ext); 
    kee = apf::DynamicMatrix(n_ext, n_ext);
    kii.zero();
    kei.zero();
    kie.zero();
    kee.zero();


    // This doesn't work. We need to get the global node number from the mesh vertices
    // the interior and exterior vectors. Still not sure how to do that.

    // Develop map from pairing of mesh nodes into the submatrix entries 
    for(i = 0; i < n_int; i++){
      for(j = 0; j < n_int; j++){
        vert_i = global_node_number(interior[i]);
        vert_j = global_node_number(interior[j]);
        k_map[std::pair<int,int>(vert_i,vert_j)] = kii.data() + n_int * i + j;
      }
    }
    for(i = 0; i < n_ext; i++){
      for(j = 0; j < n_ext; j++){
        vert_i = global_node_number(exterior[i]);
        vert_j = global_node_number(exterior[j]);
        k_map[std::pair<int,int>(vert_i,vert_j)] = kee.data() + n_ext * i + j;
      }
    }
    for(i = 0; i_ext < n_ext; i++){
      for(j = 0; j_int < n_int; j++){
        vert_i_ext = global_node_number(exterior[i_ext]);
        vert_j_int = global_node_number(interior[j_int]);
        k_map[std::pair<int,int>(vert_i_ext,vert_j_int)] = kie.data() + n_int * i_ext + j_int;
        k_map[std::pair<int,int>(vert_i_ext,vert_j_int)] = kie.data() + n_ext * j_int + i_ext;
      }
    }


    // Construct integrators for each element, assigning materials and models by to elements by model regions

    RVEvolume = 
    std::map<int, apf::Matrix3x3> mappa; // used to place kappa into IPfield by region
    auto * gmodel = mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, 3);
    while ((gent = gmi_next(gmodel, it)))
    {
      RVEvolume += gent.measure(); // Convenient time to compute
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
      // TODO: This needs to be unique_ptr-ized, as it currently isn't deleted anywhere.
      apf::Matrix3x3 * kappa_r = new apf::Matrix3x3(
        k  , 0.0, 0.0,
        0.0, k  , 0.0,
        0.0, 0.0, k  );
      mappa[tag] = *kappa_r;

      std::cout << "continuum model type: " << continuum_model->GetType()
                << "\n";
      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<EffectiveKappaIntegrator>(
              apf_primary_field, apf_primary_numbering, kappa_r, k_map);
    }
    gmi_end(gmodel, it);

    amsi::applyUniqueRegionTags(apf_mesh);
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(3);
    while((ent = apf_mesh->iterate(mesh_it)))
    {
      int tag = apf_mesh->getModelTag(apf_mesh->toModel(ent));
      apf::setMatrix(kappa, ent, 0, mappa.at(tag));
    }
  }

  EffectiveKappaStep::~EffectiveKappaStep()
  {
    apf::destroyField(apf_primary_field);
    apf::destroyField(flux);
    apf::destroyField(kappa);
    apf::destroyField(sources);
  }

  //Assemble kii, kie, kei, and kee
  //We don't need an elemental assembly, we need a vertex-pairing assembly
  void EffectiveKappaStep::Assemble()
  {
    apf::Field *coordinates = apf_mesh->getCoordinateField();
    apf::MeshIterator * mesh_region_iter = apf_mesh->begin(3);
    while (apf::MeshEntity * mesh_entity = apf_mesh->iterate(mesh_region_iter))
    {
      if (!apf_mesh->isOwned(mesh_entity))
      {
        continue;
      }
      apf::MeshElement * mlm = apf::createMeshElement(coordinates, mesh_entity);
      auto * sys = getIntegrator(mesh_entity, 0);
      //Processing this integrator directly accumulates into stiffness ubmatrices
      sys->process(mlm);  

      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(mesh_region_iter);
  }

  amsi::ElementalSystem * EffectiveKappaStep::getIntegrator(
    apf::MeshEntity * mesh_entity, int integration_point)
  {
    return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }

  void EffectiveKappaStep::FormReducedStiffnessMatrix()
  {
      kii_inverse = invert(kii);
      // k_star is currently numbered with respect to indices into 'exterior'
      k_star = kee - kei * kii_inverse * kie;
  }

  void EffectiveKappaStep::AssembleEffectiveKappa(void)
  {
    apf::Field *coordinates = apf_mesh->getCoordinateField();
    int n_ext = exterior.size();
    for(i = 0; i < n_ext; i++){
      auto vert_ext = apf_mesh.getVertex(exterior[i_ext]);
      r0 += vert_ext.evalField(coordinates);
    }
    r0 /= float(n_ext);

    apf::MeshIterator * mesh_node_iter = apf_mesh->begin(0);
    while (apf::MeshEntity * mesh_entity = apf_mesh->iterate(mesh_node_iter))
    {
    }
    for(A = 0; A < n_ext; i++){
      for(B = 0; B < n_ext; j++)
      {
        auto vert_A = apf_mesh.getVertex(exterior[A]);
        auto vert_B = apf_mesh.getVertex(exterior[B]);
        Xs_A = vert_A.evalField(coordinates) - r0;
        Xs_B = vert_B.evalField(coordinates) - r0;
        for(i = 0; i < 3; i++){
          for(j = 0; j < 3; j++){
            kappaM(i,j) += k_star[A][B] * Xs_A[i] * Xs_B[j];
          }
        }
      }
    }

  }


}