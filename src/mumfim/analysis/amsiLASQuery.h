#ifndef AMSI_LAS_QUERY_H_
#define AMSI_LAS_QUERY_H_
#include "amsiLAS.h"
#include <amsiOperators.h>
#include <cassert>
namespace amsi
{
  typedef void(LAS::*LASQueryFunc)(double&);
  struct LASNormQuery : public to_R1
  {
    LASNormQuery(LAS * l, LASQueryFunc _f)
      : las(l)
      , f(_f)
    { }
    virtual double operator()()
    {
      assert(las);
      double nrm = 0.0;
      (las->*f)(nrm);
      return nrm;
    }
  private:
    LAS * las;
    LASQueryFunc f;
  };
}
#endif
