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

extern "C" {
  extern void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  extern void dgetri_(int *N, double *A, int * LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
  extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                    double *alpha, double *A, int *lda, double *B,
                    int *ldb, double *beta, double *c, int *ldc);
}


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

    model_volume = 0.0;
    apf::MeshEntity *ent;
    auto *mesh_it = mesh->begin(3);
    while((ent = mesh->iterate(mesh_it)))
    {
      int tag = mesh->getModelTag(mesh->toModel(ent));
      apf::setMatrix(kappa, ent, 0, mappa.at(tag));
      model_volume += apf::measure(apf_mesh, ent);
    }
    std::cout << "Model volume: " << model_volume << "\n";
  
  }

  EffectiveKappaEvaluator::~EffectiveKappaEvaluator()
  {
    apf::destroyField(apf_primary_field);
    apf::destroyField(kappa);
  }

  // Populate the interior and exterior member vectors with (apf::MeshEntity*) vertices
  void EffectiveKappaEvaluator::locateVertices(void)
  {

    // Populate vert2subvert mapping and onExterior array.
    // vert2subvert[local vertex numbering] --> submatrix vertex numbering
    // Distinguish between destination numberings with onExterior.
    int number_mesh_vertices = apf_mesh->count(0);
    onExterior.reserve(number_mesh_vertices);
    onExterior.assign(number_mesh_vertices, false);
    vert2subvert.reserve(number_mesh_vertices);
    vert2subvert.assign(number_mesh_vertices,-1);
    n_int = 0, n_ext = 0;
    centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
    apf::Vector3 p;
    apf::MeshEntity *boundaryVerts[3];
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(2);
    while(ent = apf_mesh->iterate(mesh_it))
    {
      bool onBoundary = apf_mesh->countUpward(ent) == 1;
      int n_boundaryVerts = apf_mesh->getDownward(ent, 0, boundaryVerts);
      for(int ibv = 0; ibv < n_boundaryVerts; ibv++)
      {
        int ivert_number = apf::getNumber(apf_primary_numbering, 
                                          boundaryVerts[ibv], 0, 0);
        if(onBoundary)
        {
          if(vert2subvert[ivert_number] < 0) // A new, unvisited vertex
          {
            apf::getVector(coordinates, boundaryVerts[ibv], 0, p);
            centroid += p;
            onExterior[ivert_number] = true;
            vert2subvert[ivert_number] = n_ext++;
          }

        } else {
          if(vert2subvert[ivert_number] < 0) // A new, unvisited vertex
          {
            vert2subvert[ivert_number] = n_int++;
          }
        }
      }
      
    }
    centroid =  centroid / float(n_ext);
    assert(n_int + n_ext == number_mesh_vertices);

  }

  void EffectiveKappaEvaluator::Assemble(amsi::LAS *las)
  {
    locateVertices();
    //AssembleIntegratorIntoMat();
    AssembleIntegratorIntoMat_LA();
    //Solve();
    Solve_LA();
  }

  //---------------------LAPACK VERSION -----------------
  void EffectiveKappaEvaluator::Solve_LA()
  {
    assert(n_ext > 0); // Don't call this before Assemble

    /*
    double *Kii_LA_copy = new double[n_int * n_int];
    for (int i = 0; i < n_int * n_int; i++){ Kii_LA_copy[i] = Kii_LA[i];}
    double *I = new double[n_int * n_int];
    */

    double *k_star;
    if(n_int == 0){
      k_star = Kee_LA;

      printf(" k_star %dx%d :\n", n_ext, n_ext);
      for(int i=0; i<n_ext; i++){
        for(int j=0; j<n_ext; j++){
          printf("%7.4f ",k_star[i * n_ext + j]);
        }
        printf("\n");
      }
    } else {

      int *ipiv = new int[n_int];
      int info;
      int lwork = n_int * 64;
      double *work = new double[lwork];
      dgetrf_(&n_int, &n_int, Kii_LA, &n_int, ipiv, &info);
      assert(info == 0);
      dgetri_(&n_int, Kii_LA, &n_int, ipiv, work, &lwork, &info);
      assert(info == 0);

      double *KiiInverse = Kii_LA; // Kii_LA is overwritten and should be thought of as its inverse

      double alpha = 1.0, beta = 0.0;
      double *intermediate = new double[n_int * n_ext];
      char * EN = "N";
      char * TEE = "T";
      // intermediate = Inv(k^{ff}) * k^{fp}, part of second term in equation (35)
      dgemm_(EN, EN, &n_int, &n_ext, &n_int, &alpha, KiiInverse, &n_int, Kie_LA, &n_int, &beta, intermediate, &n_int);

      alpha = -1.0; beta = +1.0;
      // k* =  k^{pp} - k^{pf} * intermediate, rest of equation (35)
      dgemm_(EN, EN, &n_ext, &n_ext, &n_int, &alpha, Kei_LA, &n_ext, intermediate, &n_int, &beta, Kee_LA, &n_ext );

      k_star = Kee_LA;  // Kee is overwritten and should be thought of as k_star now.
    }

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
        //printf("%3d    %8.4f   %8.4f   %8.4f \n", ext_vert, Xsc[ext_vert][0], Xsc[ext_vert][1], Xsc[ext_vert][2]);
        ext_vert++;
      }
    }

    Km = Km * 0.0;
    double k_star_AB;
    for(int A = 0; A < n_ext; A++){
      for(int B = 0; B < n_ext; B++){
        k_star_AB = *(k_star + A * n_ext + B);
        for(int i = 0; i < 3; i++)
          //Km[i][i] += k_star_AB * Xsc[B][i] * Xsc[A][i];
          for(int j = 0; j < 3; j++)
            Km[i][j] += k_star_AB * Xsc[B][j] * Xsc[A][i];
      }
    }
    Km = Km / model_volume;

    for(int i = 0; i < 3; i++)
      printf("%8.4f  %8.4f  %8.4f\n", Km[i][0], Km[i][1], Km[i][2]);
    //std::cout << "Km: \n" <<  Km << "\n"; 
  }


  void EffectiveKappaEvaluator::AssembleIntegratorIntoMat_LA()
  {

    Kii_LA = new double[n_int*n_int];
    Kie_LA = new double[n_int*n_ext];
    Kei_LA = new double[n_ext*n_int];
    Kee_LA = new double[n_ext*n_ext];
    for(int i=0; i< n_int*n_int; i++) Kii_LA[i] = 0.0;
    for(int i=0; i< n_int*n_ext; i++) Kie_LA[i] = 0.0;
    for(int i=0; i< n_ext*n_int; i++) Kei_LA[i] = 0.0;
    for(int i=0; i< n_ext*n_ext; i++) Kee_LA[i] = 0.0;

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

      AssembleDOFs_LA(sys->getFieldNumbers(), sys->getKe());
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(mesh_region_iter);
  }


void EffectiveKappaEvaluator::AssembleDOFs_LA(const std::vector<int> & dof_numbers,
                         apf::DynamicMatrix & Ke)
  {
    
    for(int i = 0; i < 4; i++)
    {
      int i_v = dof_numbers[i];
      int i_sv = vert2subvert[i_v];
      for(int j = 0; j < 4; j++)
      {
        int j_v = dof_numbers[j];
        int j_sv = vert2subvert[j_v];
        double value = Ke(i,j);
        // These are going to need to be row major for the fortran solves
        if(onExterior[i_v]) {
          if(onExterior[j_v]) {
            *(Kee_LA + i_sv * n_ext + j_sv) += value;
          } else {
            *(Kei_LA + i_sv * n_ext + j_sv) += value;
          }
        } else {
          if(onExterior[j_v]) {
            *(Kie_LA + i_sv * n_int + j_sv) += value;
          } else {
            *(Kii_LA + i_sv * n_int + j_sv) += value;
          }
        }
      }
    }
  }


//-------------------------- PETSC VERSION -------------------------

  void EffectiveKappaEvaluator::Solve()
  {
    assert(n_ext > 0); // Don't call this before Assemble
    Mat I, KiiInverse, Product;
    MumfimPetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD, n_int, n_int, n_int, n_int, 1.0, &I));
    MumfimPetscCall(MatCreateSeqDense(PETSC_COMM_WORLD, n_int, n_int, NULL, &KiiInverse));

    // This hurts me.
    MumfimPetscCall(MatLUFactor(Kii, nullptr, nullptr, nullptr));
    MumfimPetscCall(MatMatSolve(Kii, I, KiiInverse));

    MumfimPetscCall(MatProductCreate(Kei, KiiInverse, Kie, &Product));
    MumfimPetscCall(MatProductSetType(Product, MATPRODUCT_ABC));
    MumfimPetscCall(MatProductNumeric(Product));

    MumfimPetscCall(MatAXPY(Kee, -1.0, Product, UNKNOWN_NONZERO_PATTERN));

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

// ----------------- UTILITY FUNCTIONS -----------------
  amsi::ElementalSystem * EffectiveKappaEvaluator::getIntegrator(
      apf::MeshEntity * mesh_entity,
      int /*unused integration_point */)
  {
    return constitutives[apf_mesh->toModel(mesh_entity)].get();
  }

}