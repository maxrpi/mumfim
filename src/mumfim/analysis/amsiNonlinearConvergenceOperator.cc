#include "amsiNonlinearConvergenceOperator.h"
#include <amsiOperators.h>
#include <string>
#include "amsiLASQuery.h"
#include <model_traits/CategoryNode.h>
namespace amsi {
  using mt::GetCategoryModelTraitByType;
  using mt::IModelTrait;
  using mt::IntMT;
  using mt::MTCast;
  using mt::ScalarFunctionMT;
  using mt::ScalarMT;
  static LASQueryFunc cvg_ops[] = {&LAS::GetVectorNorm, &LAS::GetSolutionNorm,
                                   &LAS::GetDotNorm};
  static LASQueryFunc ref_ops[][3] = {
      {&LAS::GetAccumVectorNorm, &LAS::GetAccumVectorNorm,
       &LAS::GetPrevVectorNorm},
      {&LAS::GetAccumSolutionNorm, &LAS::GetAccumSolutionNorm,
       &LAS::GetPrevSolutionNorm},
      {&LAS::GetAccumDotNorm, &LAS::GetAccumDotNorm, &LAS::GetPrevDotNorm}};
  to_R1* getConvergenceValueOp(int cvg_tp, LAS* las)
  {
    return new LASNormQuery(las, cvg_ops[cvg_tp]);
  }
  to_R1* getReferenceValueOp(int ref_tp, int cvg_tp, LAS* las)
  {
    if (ref_tp == 0)  // absolute
      std::cerr << "ERROR: initial reference is unimplemented." << std::endl;
    return new LASNormQuery(las, ref_ops[cvg_tp][ref_tp]);
  }
  std::unique_ptr<Convergence> createConvergenceOperator(
      const mt::CategoryNode* nd, MultiIteration* it, LAS* las)
  {
    if(nd == nullptr) {
      return nullptr;
    }
    if (nd->GetType() == "linear convergence") {
      return std::make_unique<LinearConvergence>();
    }
    if (nd->GetType() == "nonlinear iteration") {
      auto* cvg_tp =
          mt::GetCategoryModelTraitByType<IntMT>(nd, "convergence type");
      auto* ref_tp =
          mt::GetCategoryModelTraitByType<IntMT>(nd, "reference value");
      auto* eps_att = mt::GetCategoryModelTraitByType(nd, "epsilon");
      if (cvg_tp == nullptr || eps_att == nullptr || ref_tp == nullptr) {
        std::cerr << "invalid parameters provided for nonlinear iteration";
        exit(1);
      }
      // cap_att is optional
      auto* cap_att =
          mt::GetCategoryModelTraitByType<IntMT>(nd, "iteration cap");
      if (cap_att) {
        // lifetime of this iteration is linked to the lifetime of the
        // associated MultiIteration
        Iteration* stop_at_max_iters = new StopAtMaxIters((*cap_att)());
        it->addIteration(stop_at_max_iters);
      }
      to_R1* cvg_vl = getConvergenceValueOp((*cvg_tp)(), las);
      to_R1* ref_vl = getReferenceValueOp((*ref_tp)(), (*cvg_tp)(), las);
      auto* scalar_eps = MTCast<ScalarMT>(eps_att);
      if (scalar_eps) {
        auto* eps_vl = new MTConstantEpsilon((*scalar_eps)());
        return createUpdatingConvergence(it, cvg_vl, eps_vl, ref_vl);
      }
      auto* function_eps = MTCast<ScalarFunctionMT<1>>(eps_att);
      if (function_eps) {
        auto* eps_vl = new MTUpdatingEpsilon(*function_eps);
        return createUpdatingConvergence(it, cvg_vl, eps_vl, ref_vl);
      }
      std::cerr << "epsilon must be either a scalar, or scalar function with "
                   "one scalar input";
      exit(1);
    }
    return nullptr;
  }
}  // namespace amsi
