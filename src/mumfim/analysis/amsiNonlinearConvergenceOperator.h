#ifndef SIM_NONLINEAR_ANALYSIS_H_
#define SIM_NONLINEAR_ANALYSIS_H_
#include "amsiNonlinearAnalysis.h"
#include "amsiOperators.h"
#include <apf.h>
#include <model_traits/AssociatedModelTraits.h>
namespace amsi
{
  class LAS;
  std::unique_ptr<Convergence> createConvergenceOperator(
      const mt::CategoryNode *nd,
                                            MultiIteration * it,
                                            LAS * las);
  class MTUpdatingEpsilon : public R1_to_R1
  {
    protected:
      mt::ScalarFunctionMT<1> eps;
    public:
    MTUpdatingEpsilon(mt::ScalarFunctionMT<1> e)
        : eps(std::move(e))
    { }
    double operator()(double t)
    {
      return eps(t);
    }
  };
  class MTConstantEpsilon : public R1_to_R1
  {
    double eps;
    public:
    MTConstantEpsilon(double e) : eps(e) {}
    double operator()(double) {
      return eps;
    }
  };
}
#endif
