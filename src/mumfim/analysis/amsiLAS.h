#ifndef AMSI_LAS_H_
#define AMSI_LAS_H_
#include <ostream>
#include "amsiNonlinearAnalysis.h"
namespace amsi
{
  class LAS : public PerStep , public PerIter
  {
  public:
    virtual void iter() {};
    virtual void step() {};
    /// Initialize the solver and allocate memory for the local matrix portion
    virtual void Reinitialize(int,int,int,int*) {};
    virtual void Reinitialize(int,int,int) {};
    /// add an element val at given row and column in the system
    virtual void AddToMatrix(int row, int col, double value ) = 0;
    /// add a matrix of values to specific locations in the global system, useful for adding elemental system's contributions
    virtual void AddToMatrix(int num_row, int * rows,
                             int num_col, int * cols,
                             double * values) = 0;
    /// add an element in the right hand side
    virtual void AddToVector(int row, double value) = 0;
    /// add multiple elements to the right hand side
    virtual void AddToVector(int num_rows, int * rows, double * values) = 0;
    /// solve the system
    virtual void solve () = 0;
    /// zero the system
    virtual bool Zero() {return ZeroVector() && ZeroMatrix();}
    /// retrieve the solution
    virtual void GetSolution(double *&) = 0;
    virtual void GetSolutionNorm(double &) = 0;
    virtual void GetAccumSolution(double *&) = 0;
    virtual void GetAccumSolutionNorm(double &) = 0;
    virtual void GetPrevSolution(double *&) = 0;
    virtual void GetPrevSolutionNorm(double &) = 0;
    virtual bool ZeroMatrix() {return false;}
    virtual bool ZeroVector() {return false;}
    virtual void GetVector(double *&) = 0;
    virtual void SetVector(const double *) = 0;
    virtual void GetVectorNorm(double &) = 0;
    virtual void GetAccumVector(double *&) = 0;
    virtual void GetAccumVectorNorm(double &) = 0;
    virtual void GetPrevVector(double *&) = 0;
    virtual void GetPrevVectorNorm(double &) = 0;
    virtual void GetDotNorm(double &) = 0;
    virtual void GetPrevDotNorm(double &) = 0;
    virtual void GetAccumDotNorm(double &) = 0;
    virtual void PrintMatrix(std::ostream &) {};
    virtual void PrintVector(std::ostream &) {};
    virtual void PrintSolution(std::ostream &) {};
    virtual double MatrixMax() { return 0.0; }
    virtual int GlobalDOFs() { return 0; }
    virtual int LocalDOFs() { return 0; }
    virtual int LocalOffset() { return 0; }
    virtual ~LAS() {};
  };
  /*
    class VectorAccess
    {
    private:
    LAS * las;
    int rw;
    public:
    VectorAdder(LAS * l)
    : las(l)
    , rw(-1)
    { }
    VectorAdder & operator[](int r)
    {
    rw = r;
    }
    double operator+=(double val)
    {
    las->addToVector(rw,val);
    }
    };
    class MatrixAccess
    {

    };
  */
}
#endif
