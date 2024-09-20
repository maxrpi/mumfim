#include "EffectiveKappaEvaluator.h"

#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfLabelRegions.h>
#include <apfMesh.h>
#include <apf.h>
#include "mumfim/macroscale/PetscSNES.h"

#include <array>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ThermalStiffnessIntegrator.h"
#include "amsiFEA.h"
#include "gmi.h"

namespace mumfim
{
  EffectiveKappaEvaluator::EffectiveKappaEvaluator(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm comm_)
      : amsi::FEAStep(mesh, analysis_case, {}, {}, "macro", comm_)
      , constitutives()
  {
    // This makes getting vertices easier.
    apf_primary_field = apf::createLagrangeField(apf_mesh, "dummy", apf::SCALAR, 1);
    apf::zeroField(apf_primary_field);
    
    // OK, this numbering needs to be different. Or maybe we should just have
    // a reindexing array/map for each of the submatrices.
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    apf::naiveOrder(apf_primary_numbering);
    // These are convenient to have at various points.
    coordinates = apf_mesh->getCoordinateField(); 

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
          std::make_unique<ThermalStiffnessIntegrator>(
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
  
  }

  EffectiveKappaEvaluator::~EffectiveKappaEvaluator()
  {
    apf::destroyField(apf_primary_field);
    apf::destroyField(kappa);
  }

  // Populate the interior and exterior member vectors with (apf::MeshEntity*) vertices
  void EffectiveKappaEvaluator::locateVertices(void)
  {
     // Develop a vector of model faces on the boundary
    std::vector<apf::ModelEntity *> boundary_model_faces;
    auto * gmodel = apf_mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, 2); // Loop over model faces
    while ((gent = gmi_next(gmodel, it)))
    {
      gmi_set *adjacent_regions = gmi_adjacent(gmodel, gent, 3);
      if(adjacent_regions->n == 1){
        boundary_model_faces.push_back(reinterpret_cast<apf::ModelEntity *>(gent));
      }
    }
    gmi_end(gmodel, it);

    // Populate vert2subvert mapping and onExterior array.
    // vert2subvert[local vertex numbering] --> submatrix vertex numbering
    // Distinguish between destination numberings with onExterior.
    int number_mesh_vertices = apf_mesh->count(0);
    onExterior.reserve(number_mesh_vertices);
    onExterior.assign(number_mesh_vertices, false);
    vert2subvert.reserve(number_mesh_vertices);
    vert2subvert.assign(number_mesh_vertices,std::numeric_limits<int>::quiet_NaN());
    n_int = 0, n_ext = 0;
    apf::Vector3 p;
    apf::MeshEntity *boundaryVerts[3];
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(2);
    while(ent = apf_mesh->iterate(mesh_it))
    {
      apf::ModelEntity* modelEntityClassifier = apf_mesh->toModel(ent);
      int ClassifierDimension = apf_mesh->getModelType(modelEntityClassifier);
      int faceModelTag = (ClassifierDimension != 2) ?
         apf_mesh->getModelTag(modelEntityClassifier) :
         -1;
      for(auto & bmf : boundary_model_faces)
      {
        if(apf_mesh->getModelTag(bmf) == faceModelTag)
        {
          int n_boundaryVerts = apf_mesh->getDownward(ent, 0, boundaryVerts);
          for(int ibv = 0; ibv < n_boundaryVerts; ibv++)
          {
            assert(apf_mesh->getType(boundaryVerts[ibv]) == 0);
            int ivert_number = apf::getNumber(apf_primary_numbering, 
                                              boundaryVerts[ibv], 0, 0);
            if(!onExterior[ivert_number])
            {
              apf::getVector(coordinates, boundaryVerts[ibv], 0, p);
              centroid += p;
              onExterior[ivert_number] = true;
              vert2subvert[ivert_number] = n_ext++;
            }
          }
          break;
        }
      }
      
    }
    centroid =  centroid / float(n_ext);
    n_int = number_mesh_vertices - n_ext;

  }

  std::vector<int> meshVertsFromFace(apf::MeshEntity *ent)
  {

  }
  void EffectiveKappaEvaluator::Assemble(amsi::LAS *las)
  {
    locateVertices();
    AssembleIntegratorIntoMat();
    Solve();
  }

  void EffectiveKappaEvaluator::Solve()
  {
    assert(n_ext > 0); // Don't call this before Assemble
    Mat I, KiiInverse, *Product;
    MumfimPetscCall(MatCreateConstantDiagonal(analysis_comm, n_int, n_int, n_int, n_int, 1.0, &I));
    MumfimPetscCall(MatSetSizes(KiiInverse, n_int, n_int, n_int, n_int));

    // This hurts me.
    MumfimPetscCall(MatLUFactor(Kii, nullptr, nullptr, nullptr));
    MumfimPetscCall(MatMatSolve(Kii, I, KiiInverse));

    MumfimPetscCall(MatProductCreate(Kei, KiiInverse, Kie, Product));
    MumfimPetscCall(MatProductSetType(*Product, MATPRODUCT_ABC));
    MumfimPetscCall(MatProductNumeric(*Product));

    MumfimPetscCall(MatAXPY(Kee, -1.0, *Product, UNKNOWN_NONZERO_PATTERN));

    Mat &k_star = Kee;  // Kee is overwritten and should be thought of as k_star now.

    std::vector<std::array<double,3>> Xsc;
    Xsc.reserve(n_ext);
    apf::Vector3 p;
    int ext_vert = 0;
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(0);
    while((ent = apf_mesh->iterate(mesh_it)))
    {
      int vertNumber = apf::getNumber(apf_primary_numbering, ent, 0, 0);
      if(onExterior[vertNumber])
      {
        apf::getVector(coordinates, ent, 0, p);
        Xsc[ext_vert][0] = p[0] - centroid[0];
        Xsc[ext_vert][1] = p[1] - centroid[1];
        Xsc[ext_vert][2] = p[2] - centroid[2];
        ext_vert++;
      }
    }

    Km = Km * 0.0;
    double k_star_AB;
    for(int A = 0; A < n_ext; A++){
      for(int B = 0; B < n_ext; B++){
        MumfimPetscCall(MatGetValue(k_star, A, B, &k_star_AB));
        for(int i = 0; i < 3; i++)
          for(int j = 0; j < 3; j++)
            Km[i][j] += k_star_AB * Xsc[A][i] * Xsc[B][j];
      }
    }

    std::cout << "Km: \n" <<  Km << "\n"; 
  }


  amsi::ElementalSystem * EffectiveKappaEvaluator::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
    return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }

  void EffectiveKappaEvaluator::AssembleIntegratorIntoMat()
  {

    MumfimPetscCall(MatCreate(analysis_comm, &Kii));
    MumfimPetscCall(MatCreate(analysis_comm, &Kie));
    MumfimPetscCall(MatCreate(analysis_comm, &Kei));
    MumfimPetscCall(MatCreate(analysis_comm, &Kee));
    MumfimPetscCall(MatSetSizes(Kii, PETSC_DECIDE, PETSC_DECIDE, n_int, n_int));
    MumfimPetscCall(MatSetSizes(Kie, PETSC_DECIDE, PETSC_DECIDE, n_int, n_ext));
    MumfimPetscCall(MatSetSizes(Kei, PETSC_DECIDE, PETSC_DECIDE, n_ext, n_int));
    MumfimPetscCall(MatSetSizes(Kee, PETSC_DECIDE, PETSC_DECIDE, n_ext, n_ext));

    apf::MeshIterator * mesh_region_iter = apf_mesh->begin(analysis_dim);
    while (apf::MeshEntity * mesh_entity = apf_mesh->iterate(mesh_region_iter))
    {
      if (!apf_mesh->isOwned(mesh_entity))
      {
        continue;
      }
      apf::MeshElement * mlm = apf::createMeshElement(coordinates, mesh_entity);
      auto * sys = getIntegrator(mesh_entity, 0);
      sys->process(mlm);

      AssembleDOFs(sys->getFieldNumbers(), sys->getKe());
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(mesh_region_iter);
  }


void EffectiveKappaEvaluator::AssembleDOFs(const std::vector<int> & dof_numbers,
                         apf::DynamicMatrix & Ke)
  {
    
    for(int i = 0; i < 4; i++)
    {
      MumfimPetscCall(MatAssemblyBegin(Kii, MAT_FLUSH_ASSEMBLY));
      MumfimPetscCall(MatAssemblyBegin(Kie, MAT_FLUSH_ASSEMBLY));
      MumfimPetscCall(MatAssemblyBegin(Kei, MAT_FLUSH_ASSEMBLY));
      MumfimPetscCall(MatAssemblyBegin(Kee, MAT_FLUSH_ASSEMBLY));
      int i_v = dof_numbers[i];
      PetscInt i_sv = vert2subvert[i_v];
      for(int j = 0; j < 4; j++)
      {
        int j_v = dof_numbers[j];
        PetscInt j_sv = vert2subvert[j_v];
        PetscScalar value = Ke(i,j);
        if(onExterior[i_v]) {
          if(onExterior[j_v]) {
            MumfimPetscCall(MatSetValue(Kee, i_sv, j_sv, value, ADD_VALUES));
          } else {
            MumfimPetscCall(MatSetValue(Kei, i_sv, j_sv, value, ADD_VALUES));
          }
        } else {
          if(onExterior[j_v]) {
            MumfimPetscCall(MatSetValue(Kie, i_sv, j_sv, value, ADD_VALUES));
          } else {
            MumfimPetscCall(MatSetValue(Kii, i_sv, j_sv, value, ADD_VALUES));
          }
        }
      }
    }
    MumfimPetscCall(MatAssemblyEnd(Kii, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kie, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kei, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kee, MAT_FINAL_ASSEMBLY));
  }
}