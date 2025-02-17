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
#include <lapack.h>

#include "ThermalStiffnessIntegrator.h"
#include "amsiFEA.h"
#include "gmi.h"

extern "C" {
/*
  extern void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  extern void dgetri_(int *N, double *A, int * LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
  */
  extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                    double *alpha, double *A, int *lda, double *B,
                    int *ldb, double *beta, double *c, int *ldc);
}


namespace mumfim
{

  EffectiveKappaEvaluator::EffectiveKappaEvaluator(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   std::string ktf,
                                   MPI_Comm comm_)
      : amsi::FEAStep(mesh, analysis_case, {}, {}, "macro", comm_)
      , constitutives()
  {
    // Use to assign centroid. Set to true to map to origin centered unit cube.
    mapToUnitCube(apf_mesh, centroid, false);

    //apf::writeVtkFiles("mumfim_mesh", apf_mesh);
    auto tappa = read_kappa_map(ktf);

    // This makes getting vertices easier. There's probably an easier way than creating an entire field.
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
    
    // Build kappa tensors and and attach them to geometry
    auto * gmodel = mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, 3);
    while ((gent = gmi_next(gmodel, it)))
    {
      int tag = gmi_tag(gmodel, gent);

      double k = tappa[tag];
      
      apf::Matrix3x3 * kappa_r = new apf::Matrix3x3(
        k  , 0.0, 0.0,
        0.0, k  , 0.0,
        0.0, 0.0, k  );

      // Save for the kappa IPField assignment
      mappa[tag] = *kappa_r;

      constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
          std::make_unique<ThermalStiffnessIntegrator>(
              apf_primary_field, apf_primary_numbering, kappa_r);
    }
    gmi_end(gmodel, it);

    // Iterate over the mesh elements, populating a cell-based field for kappa
    // mostly for visualization later
    model_volume = 0.0;
    apf::MeshEntity *ent;
    auto *mesh_it = mesh->begin(3);
    int rcount = 0;
    while((ent = mesh->iterate(mesh_it)))
    {
      int tag = mesh->getModelTag(mesh->toModel(ent));
      apf::setMatrix(kappa, ent, 0, mappa.at(tag));
      model_volume += apf::measure(apf_mesh, ent);
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

    // Populate vert2subvert mapping and onExterior array.
    // vert2subvert[local vertex numbering] --> submatrix vertex numbering
    // Distinguish between destination numberings with onExterior.
    int number_mesh_vertices = apf_mesh->count(0);
    onExterior.assign(number_mesh_vertices, false);
    vert2subvert.assign(number_mesh_vertices,-1);

    n_int = 0, n_ext = 0;
    apf::Vector3 p;
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(2);
    while(ent = apf_mesh->iterate(mesh_it))
    {
      apf::Downward boundaryVerts;
      bool onBoundary = apf_mesh->countUpward(ent) == 1;
      int n_boundaryVerts = apf_mesh->getDownward(ent, 0, boundaryVerts);

      if(onBoundary)
      {
        for(int ibv = 0; ibv < n_boundaryVerts; ibv++)
        {
          int ivert_number = apf::getNumber(apf_primary_numbering, 
                                            boundaryVerts[ibv], 0, 0);
          if(vert2subvert[ivert_number] < 0) // A new, unvisited vertex
          {
            onExterior[ivert_number] = true;
            vert2subvert[ivert_number] = n_ext++;
          }
        }
      }
    }
    apf_mesh->end(mesh_it);

    mesh_it = apf_mesh->begin(2);
    while(ent = apf_mesh->iterate(mesh_it))
    {
      apf::Downward boundaryVerts;
      bool onBoundary = apf_mesh->countUpward(ent) == 1;
      int n_boundaryVerts = apf_mesh->getDownward(ent, 0, boundaryVerts);
      if(!onBoundary)
      {
        for(int ibv = 0; ibv < n_boundaryVerts; ibv++)
        {
          int ivert_number = apf::getNumber(apf_primary_numbering, 
                                            boundaryVerts[ibv], 0, 0);
          if(vert2subvert[ivert_number] < 0) // A new, unvisited vertex
          {
            vert2subvert[ivert_number] = n_int++;
          }
        }
      }
    }
    apf_mesh->end(mesh_it);
    assert(n_int + n_ext == number_mesh_vertices);

  }

  void EffectiveKappaEvaluator::Assemble(amsi::LAS *las)
  {
    locateVertices();
    AssembleIntegratorIntoMat_LA();
    Solve_LA();
  }

  void EffectiveKappaEvaluator::sanityCheck1(double *k_star,
                    std::vector<std::array<double,3>> &Xsc,
                    int n_ext, int n_int, bool colMajor)
  {
    double sum[3] = {0.0, 0.0, 0.0};
    for(int A=0; A < n_ext + n_int; A++)
    {
      if(!onExterior[A]) continue;
      int sA = vert2subvert[A];
      for(int B=0; B < n_ext + n_int; B++)
      {
        if(!onExterior[B]) continue;
        int sB = vert2subvert[B];
        //printf("%d, %d: %f \n", sA, sB, colMajor ? *(k_star + sA + N * sB) : *(k_star + sA * N + sB));
        for(int i=0; i < 3; i++)
        {
          sum[i] += colMajor 
            ? *(k_star + sA + n_ext * sB) * Xsc[sA][i]
            : *(k_star + sA * n_ext + sB) * Xsc[sA][i];
        }
      }
    }
    std::cout << " sum of kstar_AB * X_A:\ ";
    for(int i=0; i < 3; i++)
    {
      std::cout  << "  " << sum[i] << " \n";
    }
    std::cout << "\n";
  }

  bool isTranspose(int N, int M, double *A, double *B, double tol=1.0e-16)
  {
    for(int i=0; i<M; i++)
      for(int j=0; j<N; j++)
        if(abs(*(A + i * N + j) - *(B + j * M + i)) > tol) return false;
    return true;
  }

  //---------------------LAPACK VERSION -----------------
  void EffectiveKappaEvaluator::Solve_LA()
  {
    assert(n_ext > 0); // Don't call this before Assemble
    printf("n_int = %d, n_ext = %d\n", n_int, n_ext);

    double *k_star;
    if(n_int == 0){
      k_star = Kee_LA;

    } else {

      int *ipiv = new int[n_int];
      int info;
      dgetrf_(&n_int, &n_int, Kii_LA, &n_int, ipiv, &info);
      assert(info == 0);
      double work_query[1];
      int lwork_query = -1;
      dgetri_(&n_int, Kii_LA, &n_int, ipiv, work_query, &lwork_query, &info);
      int lwork = int(work_query[0]);
      double *work = new double[lwork];
      dgetri_(&n_int, Kii_LA, &n_int, ipiv, work, &lwork, &info);
      assert(info == 0);

      double *KiiInverse = Kii_LA; // Kii_LA is overwritten and should be thought of as its inverse

      double alpha = 1.0, beta = 0.0;
      double *intermediate = new double[n_int * n_ext];
      char EN = 'N';
      char TEE = 'T';

      // dgemm has this form:
      // dgemm_(transA, transB, M, N, K,   alpha, A_MxK, ldA_MxK,   B_KxN, ldB_KxN,    beta, C_MxN, ldC_MxN)

      // intermediate = Inv(k^{ff}) * k^{fp}, part of second term in equation (35)
      dgemm_(&EN, &EN, &n_int, &n_ext, &n_int,   &alpha, KiiInverse, &n_int,   Kie_LA, &n_int,   &beta, intermediate, &n_int);

      alpha = -1.0; beta = +1.0;
      // k* =  k^{pp} - k^{pf} * intermediate, rest of equation (35)
      dgemm_(&EN, &EN, &n_ext, &n_ext, &n_int,   &alpha, Kei_LA, &n_ext,   intermediate, &n_int,   &beta, Kee_LA, &n_ext );

      k_star = Kee_LA;  // Kee is overwritten and should be thought of as k_star now.
    }

    std::vector<std::array<double,3>> Xsc;
    Xsc.assign(n_ext, {99.99, 99.99, 99.99});
    apf::Vector3 p;
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(0);
    while((ent = apf_mesh->iterate(mesh_it)))
    {
      int vertNumber = apf::getNumber(apf_primary_numbering, ent, 0, 0);
      if(onExterior[vertNumber])
      {
        apf::getVector(coordinates, ent, 0, p);
        Xsc[vert2subvert[vertNumber]][0] = p[0] - centroid[0];
        Xsc[vert2subvert[vertNumber]][1] = p[1] - centroid[1];
        Xsc[vert2subvert[vertNumber]][2] = p[2] - centroid[2];
      }
    }
    //sanityCheck1(k_star, Xsc, n_ext, n_int, false);

    apf::Matrix3x3 Km(0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    
    double k_star_AB;
    for(int A = 0; A < n_ext + n_int; A++){
      if(!onExterior[A]) continue;
      int sA = vert2subvert[A];
      assert(sA < n_ext);
      for(int B = 0; B < n_ext + n_int; B++){
        if(!onExterior[B]) continue;
        int sB = vert2subvert[B];
        assert(sB < n_ext);
        k_star_AB = k_star[sA * n_ext + sB]; 
        for(int i = 0; i < 3; i++)
          for(int j = 0; j < 3; j++)
            Km[i][j] += k_star_AB * Xsc[sB][j] * Xsc[sA][i];
      }
    }
   
    Km = Km / model_volume;

    std::cout << "MacroKappa: ";
    printf("%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
            centroid[0], centroid[1], centroid[2],
            Km[0][0], Km[1][1], Km[2][2], Km[0][1], Km[2][1], Km[0][2]);
    //std::cout << "MacroKappa: \n";
    //for(int i = 0; i < 3; i++)
      //printf("%8.4f  %8.4f  %8.4f\n", Km[i][0], Km[i][1], Km[i][2]);
  
    apf::destroyNumbering(apf_primary_numbering);
    delete(Kii_LA);
    delete(Kie_LA);
    delete(Kei_LA);
    delete(Kee_LA);
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
    
    for(int i = 0; i < dof_numbers.size(); i++)
    {
      int i_v = dof_numbers[i];
      int i_sv = vert2subvert[i_v];
      for(int j = 0; j < dof_numbers.size(); j++)
      {
        int j_v = dof_numbers[j];
        int j_sv = vert2subvert[j_v];
        double value = Ke(i,j);
      
        if(onExterior[i_v]) {
          assert(i_sv < n_ext);
          if(onExterior[j_v]) {
            assert(j_sv < n_ext);
            *(Kee_LA + i_sv * n_ext + j_sv) += value; 
          } else {
            assert(j_sv < n_int);
            *(Kie_LA + i_sv * n_int + j_sv) += value; 
          }
        } else {
          assert(i_sv < n_int);
          if(onExterior[j_v]) {
            assert(j_sv < n_ext);
            *(Kei_LA + i_sv * n_ext + j_sv) += value;
          } else {
            assert(j_sv < n_int);
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

    Mat *k_star;
    if(n_int == 0){
      k_star = &Kee;

    } else {
      Mat I, KiiInverse, Product;
      MumfimPetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD, n_int, n_int, n_int, n_int, 1.0, &I));
      // This hurts me.
      MumfimPetscCall(MatSetType(I, MATDENSE));
      MumfimPetscCall(MatCreateSeqDense(PETSC_COMM_WORLD, n_int, n_int, NULL, &KiiInverse));
      MumfimPetscCall(MatSetType(KiiInverse, MATDENSE));

      IS isrow, iscol;
      MatFactorInfo mfinfo;
      MumfimPetscCall(MatGetOrdering(Kii, MATORDERINGNATURAL, &isrow, &iscol));
      MumfimPetscCall(MatLUFactor(Kii, isrow, iscol, &mfinfo));
      MumfimPetscCall(MatMatSolve(Kii, I, KiiInverse));

      MumfimPetscCall(MatProductCreate(Kei, KiiInverse, Kie, &Product));
      MumfimPetscCall(MatProductSetType(Product, MATPRODUCT_ABC));
      MumfimPetscCall(MatProductNumeric(Product));

      MumfimPetscCall(MatAXPY(Kee, -1.0, Product, UNKNOWN_NONZERO_PATTERN));

      k_star = &Kee;  // Kee is overwritten and should be thought of as k_star now.
    }

    std::vector<std::array<double,3>> Xsc;
    Xsc.assign(n_ext, {99.99, 99.99, 99.99});
    apf::Vector3 p;
    apf::MeshEntity *ent;
    auto *mesh_it = apf_mesh->begin(0);
    while((ent = apf_mesh->iterate(mesh_it)))
    {
      int vertNumber = apf::getNumber(apf_primary_numbering, ent, 0, 0);
      if(onExterior[vertNumber])
      {
        apf::getVector(coordinates, ent, 0, p);
        Xsc[vert2subvert[vertNumber]][0] = p[0] - centroid[0];
        Xsc[vert2subvert[vertNumber]][1] = p[1] - centroid[1];
        Xsc[vert2subvert[vertNumber]][2] = p[2] - centroid[2];
      }
    }

    apf::Matrix3x3 Km(0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);

    double k_star_AB;
    for(int A = 0; A < n_ext + n_int; A++){
      if(!onExterior[A]) continue;
      PetscInt sA = vert2subvert[A];
      for(int B = 0; B < n_ext + n_int; B++){
        if(!onExterior[B]) continue;
        PetscInt sB = vert2subvert[B];
        MumfimPetscCall(MatGetValue(*k_star, sA, sB, &k_star_AB));
        for(int i = 0; i < 3; i++)
          for(int j = 0; j < 3; j++)
            Km[i][j] += k_star_AB * Xsc[sA][i] * Xsc[sB][j];
      }
    }
    Km = Km / model_volume;

    std::cout << "MacroKappa: ";
    printf("%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
            centroid[0], centroid[1], centroid[2],
            Km[0][0], Km[1][1], Km[2][2], Km[0][1], Km[2][1], Km[0][2]);
    
    apf::destroyNumbering(apf_primary_numbering);
    delete(Kii_LA);
    delete(Kie_LA);
    delete(Kei_LA);
    delete(Kee_LA);
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
    MumfimPetscCall(MatSetFromOptions(Kii));
    MumfimPetscCall(MatSetFromOptions(Kie));
    MumfimPetscCall(MatSetFromOptions(Kei));
    MumfimPetscCall(MatSetFromOptions(Kee));


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
    MumfimPetscCall(MatAssemblyBegin(Kii, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyBegin(Kie, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyBegin(Kei, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyBegin(Kee, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kii, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kie, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kei, MAT_FINAL_ASSEMBLY));
    MumfimPetscCall(MatAssemblyEnd(Kee, MAT_FINAL_ASSEMBLY));
  }


void EffectiveKappaEvaluator::AssembleDOFs(const std::vector<int> & dof_numbers,
                         apf::DynamicMatrix & Ke)
  {
    
    for(int i = 0; i < dof_numbers.size(); i++)
    {
      int i_v = dof_numbers[i];
      PetscInt i_sv = vert2subvert[i_v];
      for(int j = 0; j < dof_numbers.size(); j++)
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

  void mapToUnitCube(apf::Mesh *mesh, double *centroid_, bool scale)
  {
    double centroid[3];
    apf::Vector3 p;
    double min[3] = {1.0e30, 1.0e30, 1.0e30}, max[3] = {-1.0e30, -1.0e30, -1.0e30};

    apf::Field *coordinates = mesh->getCoordinateField();
    //apf::Numbering *numbering = apf::createNumbering(coordinates);

    apf::MeshEntity *ent;
    auto *mesh_it = mesh->begin(0);
    while(ent = mesh->iterate(mesh_it))
    {
      apf::getVector(coordinates, ent, 0, p);
      for(int j=0; j<3; j++){
        min[j] = p[j] < min[j] ?  p[j] : min[j];
        max[j] = p[j] > max[j] ?  p[j] : max[j];
      }
    }     
    printf("MIN: %7.4f  %7.4f  %7.4f\n", min[0], min[1], min[2]);
    printf("MAX: %7.4f  %7.4f  %7.4f\n", max[0], max[1], max[2]);
    for(int j=0; j<3; j++){ centroid[j] = (max[j] + min[j]) / 2.0;}
    printf("CENTROID: %7.4f  %7.4f  %7.4f\n", centroid[0], centroid[1], centroid[2]);
    if(!scale && centroid_ != nullptr){
      centroid_[0] = centroid[0];
      centroid_[1] = centroid[1];
      centroid_[2] = centroid[2];
    }

    if(scale)
    {
      mesh->end(mesh_it);

      mesh_it = mesh->begin(0);
      while(ent = mesh->iterate(mesh_it))
      {
        apf::getVector(coordinates, ent, 0, p);
        for(int j=0; j<3; j++){
          p[j] = (p[j] - centroid[j]) / (max[j] - min[j]);
        }
        apf::setVector(coordinates, ent, 0, p);
      }
      mesh->end(mesh_it);
      // Update with new centroid.
      centroid_[0] = 0.0; 
      centroid_[1] = 0.0; 
      centroid_[2] = 0.0; 
    }

  }

  std::map<int, double> read_kappa_map(std::string kappa_tag_filename){
    std::map<int, double> tappa;
    FILE *fp = fopen(kappa_tag_filename.c_str(), "r");
    int tag_; float kappa_;
    while(fscanf(fp, "%d, %f", &tag_, &kappa_) != EOF){
      tappa[tag_] = (double) kappa_;
    }
    fclose(fp);
    return tappa;
  }

}