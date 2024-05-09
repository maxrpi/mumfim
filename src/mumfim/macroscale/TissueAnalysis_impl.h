#include <amsiNonlinearConvergenceOperator.h>
#include <model_traits/AssociatedModelTraits.h>
#include <model_traits/ModelTraits.h>
namespace mumfim
{
  template <typename O>
  void buildLASConvergenceOperators(const mt::CategoryNode * solution_strategy,
                                    amsi::MultiIteration * it,
                                    amsi::LAS * las,
                                    O out)
  {
    auto convergence_operators_cat =
        mt::GetCategoriesByType(solution_strategy, "convergence operator");
    for (const auto * convergence_operator_cat : convergence_operators_cat)
    {
      for (const auto & convergence_operator :
           convergence_operator_cat->GetCategoryNodes())
      {
        auto * cvg =
            amsi::createConvergenceOperator(&convergence_operator, it, las)
                .release();
        if (cvg != nullptr)
        {
          *out++ = cvg;
        }
      }
    }
  }
}  // namespace mumfim
