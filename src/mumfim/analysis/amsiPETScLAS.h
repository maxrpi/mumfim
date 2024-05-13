#ifndef PETSCITERATIVESOLVER_H_
#define PETSCITERATIVESOLVER_H_
#include "amsiLAS.h"
#include <petscksp.h>
#include <petscmat.h>
#include <vector>
namespace amsi
{
  class PetscLAS : public LAS
  {
  public:
    PetscLAS();
    PetscLAS(int,int);
    void iter();
    void step();
    void Reinitialize(int, int, int, int*);
    void Reinitialize(int, int, int);
    void AddToMatrix(int,int,double);
    void AddToMatrix(int, const int *,int, const int *, const double *);
    void AddToVector(int,double);
    void AddToVector(int, const int *, const double *);
    void solve();
    bool Zero();
    bool ZeroMatrix();
    bool ZeroVector();
    void GetVector(double *&);
    void SetVector(const double *);
    void GetVectorNorm(double &);
    void GetAccumVector(double *&);
    void GetAccumVectorNorm(double &);
    void GetPrevVector(double *&);
    void GetPrevVectorNorm(double &);
    void GetDotNorm(double &);
    void GetPrevDotNorm(double &);
    void GetAccumDotNorm(double &);
    void GetSolution(double *&);
    void GetSolutionNorm(double &);
    void GetAccumSolution(double *&);
    void GetAccumSolutionNorm(double &);
    void GetPrevSolution(double *&);
    void GetPrevSolutionNorm(double &);
    void PrintMatrix(std::ostream &);
    void PrintVector(std::ostream &);
    void PrintSolution(std::ostream &);
    double MatrixMax();
    int GlobalDOFs() { return globalNumEqs; }
    int LocalDOFs() { return vec_high - vec_low; }
    int LocalOffset() { return vec_low; }
    Mat GetMatrix() { return A; }
    Vec GetVector() { return b_i; }
    Vec GetSolutionVector() { return x; }
    ~PetscLAS();
  private:
    bool isVectorAssembled();
    void freeMem();
    Mat A;     // matrix
    Vec x_im;   // previous solution
    Vec x_i;     // current solution
    Vec x;
    Vec b_im;   // previous vector
    Vec b_i;     // current vector
    Vec b;
    Vec w;     // work vector
    double * x_arr;
    double * b_arr;
    double * b_i_arr;
    int globalNumEqs;
    int vec_low;
    int vec_high;
    int mat_low;
    int mat_high;
    bool b_assembled, b_addMode;
    KSP solver;
  };
}
#endif
