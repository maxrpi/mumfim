#ifndef AMSI_DUMMY_LAS_H_
#define AMSI_DUMMY_LAS_H_
#include "amsiLAS.h"
#include "MatrixUtil.h"
#include <cmath>
namespace amsi
{
  class DummyLAS : public LAS
  {
  private:
    double ** mtx;
    double * vc;
    double * sl;
    int dof;
  public:
    DummyLAS(int dofs)
      : mtx(allocate_matrix(dofs))
      , vc(new double[dofs])
      , sl(new double[dofs])
      , dof(dofs)
    {
      memset(vc,0.0,sizeof(double)*dof);
      memset(sl,0.0,sizeof(double)*dof);
    }
    ~DummyLAS()
    {
      free_matrix(mtx);
      delete [] vc;
      delete [] sl;
    }
    virtual void AddToMatrix(int rw, int cl, double vl)
    {
      mtx[rw][cl] += vl;
    }
    virtual void AddToMatrix(int rw_cnt, int * rws,
                             int cl_cnt, int * cls,
                             double * vls)
    {
      for(int ii = 0; ii < rw_cnt; ++ii)
        for(int jj = 0; jj < cl_cnt; ++jj)
          mtx[rws[ii]][cls[jj]] += vls[rw_cnt*ii + jj];
    }
    virtual void AddToVector(int rw, double vl)
    {
      vc[rw] += vl;
    }
    virtual void AddToVector(int rw_cnt, int * rws, double * vls)
    {
      for(int ii = 0; ii < rw_cnt; ++ii)
        vc[rws[ii]] += vls[ii];
    }
    virtual void solve () {};
    virtual void GetSolution(double *& sol)
    {
      sol = sl;
    }
    virtual void GetSolutionNorm(double & nrm)
    {
      nrm = 0.0;
      for(int ii = 0; ii < dof; ii++)
        nrm += sl[ii] * sl[ii];
      nrm = sqrt(nrm);
    }
    // stubb
    virtual void GetAccumSolutionNorm(double &)
    {   }
    virtual bool ZeroMatrix()
    {
      for(int ii = 0; ii < dof; ++ii)
        memset(mtx[ii],0,dof*sizeof(double));
      return true;
    }
    virtual bool ZeroVector()
    {
      memset(vc,0,dof*sizeof(double));
      return true;
    }
    virtual void GetVector(double *& vec)
    {
      vec = vc;
    }
    // stubb
    virtual void GetAccumVector(double *&)
    {    }
    virtual void SetVector(const double * vec)
    {
      memcpy(&vc[0],&vec[0],dof*sizeof(double));
    }
    virtual void GetVectorNorm(double & nrm)
    {
      nrm = 0.0;
      for(int ii = 0; ii < dof; ii++)
        nrm += vc[ii] * vc[ii];
      nrm = sqrt(nrm);
    }
    // stubb
    virtual void GetAccumVectorNorm(double & )
    {   }
    // stubb
    virtual void GetDotNorm(double &)
    {   }
  };
}
#endif
