#include "amsiPETScLAS.h"
#include <amsiMPI.h>
#include <assert.h>
#include <iostream>
#include <limits>
#include <amsiVerbosity.h>
namespace amsi
{
  void PetscLAS::Reinitialize(int num_local_unknowns,
                              int num_global_unknowns,
                              int global_offset,
                              int * )
  {
    Reinitialize(num_local_unknowns,num_global_unknowns,global_offset);
  }
  void PetscLAS::Reinitialize(int num_local_unknowns,
                              int num_global_unknowns,
                              int global_offset)
  {
    int rsz_root = num_local_unknowns != vec_high - vec_low || global_offset != vec_low;
    int rsz = 0;
    MPI_Allreduce(&rsz_root,&rsz,1,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
    if(rsz)
    {
      if(x_arr)
        freeMem();
      int N = globalNumEqs = num_global_unknowns;
      int n = num_local_unknowns;
      x_arr = new double[n];
      b_arr = new double[n];
      b_i_arr = new double[n];
      VecCreateMPI(PETSC_COMM_WORLD,n,N,&b_i);
      VecSetOption(b_i, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
      MatCreateAIJ(PETSC_COMM_WORLD,n,n,N,N,300,PETSC_NULL,300,PETSC_NULL,&A);
      //MatMPIAIJSetPreallocation(A,300,PETSC_NULL,300,PETSC_NULL); // This is done in previous line
      MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      KSPCreate(PETSC_COMM_WORLD,&solver);
      // create other vectors...
      VecDuplicate(b_i,&x_i);
      VecDuplicate(b_i,&x_im);
      VecDuplicate(b_i,&b_im);
      VecDuplicate(b_i,&b);
      VecDuplicate(b_i,&x);
      VecDuplicate(b_i,&w);
      VecGetOwnershipRange(b_i,&vec_low,&vec_high);
      MatGetOwnershipRange(A,&mat_low,&mat_high);
      AMSI_V1(std::cout << "Local equations = " << n << ", Global Equations = " << N << std::endl;)
      AMSI_V1(std::cout << "Local equations = " << n << ", Global Equations = " << N << std::endl;)
      AMSI_V1(std::cout << "Vector ownership range: " << vec_low << "-" << vec_high << std::endl;)
      AMSI_V1(std::cout << "Matrix ownership range: " << mat_low << "-" << mat_high << std::endl;)
      AMSI_V1(std::cout << "(Re)initialized the matrix, vectors and the solver" << std::endl;)
    }
  }
  void PetscLAS::iter()
  {
    // move into previous iteration vecs
    VecCopy(b_i,b_im);
    VecCopy(x_i,x_im);
  }
  void PetscLAS::step()
  {
    // no iteration values anymore
    VecZeroEntries(b_i);
    VecZeroEntries(x_i);
    VecZeroEntries(b_im);
    VecZeroEntries(x_im);
    // zero accumulated values as well, they're already reflected in the analysis
    VecZeroEntries(b);
    VecZeroEntries(x);
  }
  /**
   *@brief Add a value to the current value at the given location in the matrix in the Linear System to be solved.
   *
   *@param[in] The row at which to ADD the value (1-indexed).
   *@param[in] The column at which to ADD the value (1-indexed).
   *@param[in] The value to insert into the matrix.
   */
  void PetscLAS::AddToMatrix(int row, int col, double value)
  {
    MatSetValues(A,1,&row,1,&col,&value,ADD_VALUES);
  }
  void PetscLAS::AddToMatrix(int num_rows, int * rows, int num_cols, int * cols, double * values)
  {
    MatSetValues(A,num_rows,rows,num_cols,cols,values,ADD_VALUES);
  }
  /**
   *@brief Add a value to the current value at the given location in the vector in the Linear System to be solved.
   *
   *@param[in] The row of the RHS vector into which to ADD the value.
   *@param[in] The value to ADD to the vector.
   *
   *@todo Construct a test case.
   */
  void PetscLAS::AddToVector(int row , double value)
  {
    if(!b_addMode)
    {
      VecAssemblyBegin(b_i);
      VecAssemblyEnd(b_i);
      b_addMode = true;
    }
    VecSetValue(b_i,row,static_cast<PetscScalar>(value),ADD_VALUES);
    b_assembled = false;
  }
  void PetscLAS::AddToVector(int num_rows, int * rows, double * values)
  {
    if(!b_addMode)
    {
      VecAssemblyBegin(b_i);
      VecAssemblyEnd(b_i);
      b_addMode = true;
    }
    VecSetValues(b_i,num_rows,rows,static_cast<PetscScalar*>(values),ADD_VALUES);
    b_assembled = false;
  }
  /**
   *@brief Solve the Linear System.
   */
  void PetscLAS::solve()
  {
    if(!isVectorAssembled())
    {
      VecAssemblyBegin(b_i);
      VecAssemblyEnd(b_i);
      b_assembled = true;
    }
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    KSPSetOperators(solver,A,A);
    KSPSetFromOptions(solver);
    //PetscObjectDereference((PetscObject)A);
    //PetscObjectDereference((PetscObject)A);
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    //VecView(b_i, PETSC_VIEWER_STDOUT_WORLD);
    // solve the system
    KSPSolve(solver,b_i,x_i);
    //PetscViewer vwr;
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"delta_u_i",&vwr);
    //VecView(x_i,vwr);
    //PetscViewerDestroy(&vwr);
    // accumulate values
    VecAYPX(b,1.0,b_i);
    VecAYPX(x,1.0,x_i);
  }
  /**
   *@brief Zero out the Matrix and Vector associated with the linear system, without changing their underlying nonzero structures.
   *
   *@todo Construct a test case.
   */
  bool PetscLAS::Zero()
  {
    return (ZeroVector() && ZeroMatrix());
  }
  /**
   *@brief Retrieve a pointer to an array of doubles containing the solution to the linear system.
   *
   *@todo Construct a test case. Prevent the code from being executed prior to a solution being produced.
   */
  void PetscLAS::GetSolution(double *& sol)
  {
    PetscScalar * X;
    VecGetArray(x_i,&X);
    int n = vec_high-vec_low;
    memcpy(x_arr,(double*)X,n*sizeof(PetscScalar));
    VecRestoreArray(x_i,&X);
    sol = &x_arr[0];
  }
  /**
   *@brief Zero out the associated Matrix, should be collective on the Matrix.
   *
   *@return Whether or not the matrix was successfully zeroed.
   *
   *@todo Construct a test case in a parallel environment in order to assure the operation is collective.
   */
  bool PetscLAS::ZeroMatrix()
  {
    // logically collective on A
    return (MatZeroEntries(A) == 0);
  }
  /**
   *@brief Zero out the associated Vector, should be collective on the Vector.
   *
   *@return Whether or not the vector was successfully zeroed.
   *
   *@todo Construct a test case in a parallel environment in order to asufe the operation is collective.
   */
  bool PetscLAS::ZeroVector()
  {
    if(b_addMode)
    {
      VecAssemblyBegin(b_i);
      VecAssemblyEnd(b_i);
    }
    b_assembled = false;
    VecZeroEntries(b_i);
    return true;
  }
  /**
   *@brief Retrieve an array of doubles which is a serial copy of the associated vector.
   *
   *@return A pointer to an array of \c doubles representing the associated vector.
   *
   *@todo Construct a test case.
   */
  void PetscLAS::GetVector(double * & vc)
  {
    PetscScalar * Bi;
    VecGetArray(b_i,&Bi);
    int n = vec_high-vec_low;
    memcpy(b_i_arr,(double*)Bi,n*sizeof(PetscScalar));
    VecRestoreArray(b_i,&Bi);
    vc = &b_i_arr[0];
  }
  void PetscLAS::SetVector(const double * vec)
  {
    if(!b_addMode)
    {
      VecAssemblyBegin(b_i);
      VecAssemblyEnd(b_i);
    }
    // This should allow us to update the entire vector at once without the need to generate the array as below, still UNTESTED
    VecSetBlockSize(b_i,globalNumEqs);
    int zero = 0;
    VecSetValuesBlocked(b_i,1,&zero,vec,INSERT_VALUES);
    b_assembled = false;
    b_addMode = false;
  }
  void PetscLAS::GetVectorNorm(double & norm)
  {
    if(b_i)
    {
      if(!isVectorAssembled())
      {
        VecAssemblyBegin(b_i);
        VecAssemblyEnd(b_i);
        b_assembled = true;
      }
      VecNorm(b_i,NORM_2,&norm);
    }
    else
      norm = 0.0;
  }
  void PetscLAS::GetAccumVector(double * & vc)
  {
    PetscScalar * B;
    VecGetArray(b,&B);
    int n = vec_high-vec_low;
    memcpy(b_arr,(double*)B,n*sizeof(PetscScalar));
    VecRestoreArray(b,&B);
    vc = &b_arr[0];
  }
  void PetscLAS::GetAccumVectorNorm(double & nrm)
  {
    if(b)
      VecNorm(b,NORM_2,&nrm);
    else
      nrm = 0.0;
  }
  void PetscLAS::GetPrevVector(double *&)
  {
    // UNIMPLEMENTED
  }
  void PetscLAS::GetPrevVectorNorm(double & nrm)
  {
    if(b_im)
      VecNorm(b_im,NORM_2,&nrm);
    else
      nrm = 0.0;
  }
  void PetscLAS::GetDotNorm(double & nrm)
  {
    if(b_i && x_i)
    {
      VecDot(b_i,x_i,&nrm);
      nrm = fabs(nrm);
    }
    else
      nrm = 0.0;
  }
  void PetscLAS::GetPrevDotNorm(double & nrm)
  {
    if(b_im && x_im)
    {
      VecDot(b_im,x_im,&nrm);
      nrm = fabs(nrm);
    }
    else
      nrm = 0.0;
  }
  void PetscLAS::GetAccumDotNorm(double & nrm)
  {
    if(b && x)
    {
      VecDot(b,x,&nrm);
      nrm = fabs(nrm);
    }
    else
      nrm = 0.0;
  }
  void PetscLAS::GetSolutionNorm(double & nrm)
  {
    if(x_i)
      VecNorm(x_i,NORM_2,&nrm);
    else
      nrm = 0.0;
  }
  void PetscLAS::GetAccumSolution(double *&)
  {
    // UNIMPLEMENTED
  }
  void PetscLAS::GetAccumSolutionNorm(double & nrm)
  {
    if(x)
      VecNorm(x,NORM_2,&nrm);
    else
      nrm = 0.0;
  }
  void PetscLAS::GetPrevSolution(double *&)
  {
    // UNIMPLEMENTED
  }
  void PetscLAS::GetPrevSolutionNorm(double & nrm)
  {
    if(x_im)
      VecNorm(x_im,NORM_2,&nrm);
    else
      nrm = 0.0;
  }
  PetscLAS::~PetscLAS()
  {
    freeMem();
  }
  void PetscLAS::freeMem()
  {
    KSPDestroy(&solver);
    MatDestroy(&A);
    VecDestroy(&x);
    VecDestroy(&x_i);
    VecDestroy(&x_im);
    VecDestroy(&b);
    VecDestroy(&x);
    VecDestroy(&b_i);
    VecDestroy(&b_im);
    VecDestroy(&w);
    delete [] x_arr;
    delete [] b_arr;
    delete [] b_i_arr;
  }
  PetscLAS::PetscLAS()
    : x_arr(NULL)
    , b_arr(NULL)
    , b_i_arr(NULL)
    , globalNumEqs(0)
    , vec_low(0)
    , vec_high(0)
    , mat_low(0)
    , mat_high(0)
    , b_assembled(true)
    , b_addMode(true)
  { }
  /**
   *@brief Constructor.
   *
   *@param[in] The number of global equations across the linear system.
   *@param[in] The number of local equations associated with the local system.
   */
  PetscLAS::PetscLAS ( int N, int n )
    : A()
    , x_im()
    , x_i()
    , x()
    , b_im()
    , b_i()
    , b()
    , w()
    , x_arr(NULL)
    , b_arr(NULL)
    , b_i_arr(NULL)
    , globalNumEqs(N)
    , vec_low(0)
    , vec_high(0)
    , mat_low(0)
    , mat_high(0)
    , b_assembled(true)
    , b_addMode(true)
    , solver()
  {
    int ffst = 0;
    MPI_Scan(&n,&ffst,1,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
    ffst -= n; // scan is inclusive, remove the local dofs from the offset
    Reinitialize(n,N,ffst);
  }
  /**
   * @brief Print the associated matrix to the specified output.
   *
   * @param[in] Currently this does not effect where the matrix is printed.
   */
  void PetscLAS::PrintMatrix(std::ostream &)
  {
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  }
  /**
   *@brief Print the associated vector to standard output.
   *
   *@param[in] Currently this does not effect where the vector is printed.
   */
  void PetscLAS::PrintVector(std::ostream &)
  {
    VecAssemblyBegin(b_i);
    VecAssemblyEnd(b_i);
    VecView(b_i, PETSC_VIEWER_STDOUT_WORLD);
  }
  void PetscLAS::PrintSolution(std::ostream &)
  {
    VecAssemblyBegin(x_i);
    VecAssemblyEnd(x_i);
    VecView(x_i,PETSC_VIEWER_STDOUT_WORLD);
  }
  double PetscLAS::MatrixMax()
  {
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    double rslt = 0.0;
    PetscInt idx[globalNumEqs];
    MatGetRowMax(A,w,idx);
    VecMax(w,NULL,&rslt);
    return rslt;
  }
  bool PetscLAS::isVectorAssembled( )
  {
    MPI_Comm cm;
    PetscObjectGetComm((PetscObject)b_i,&cm);
    int assembled = b_assembled;
    // not valid MPI to do min on bool type
    return comm_min(assembled,cm);
  }
}
