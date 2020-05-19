#include "bioFiberRVEAnalysis.h"
#include "bioFiberRVEAnalysisExplicit.h"
#include "bioFiberRVEAnalysisStaticImplicit.h"
#include <memory>

namespace bio {
  std::unique_ptr<FiberRVEAnalysis> createFiberRVEAnalysis(
      std::unique_ptr<FiberNetwork> fn,
      std::unique_ptr<MicroSolutionStrategy> ss,
      std::shared_ptr<void> workspace) {

    auto analysis = std::unique_ptr<FiberRVEAnalysis>{nullptr};
    if(fn == nullptr || ss == nullptr)
    {
      std::cerr<<"You must pass in a valid fiber network and solution strategy poiter!"<<std::endl;
      std::abort();
    }
    else if (ss->slvrType == SolverType::Explicit) {
       analysis.reset(new FiberRVEAnalysisExplicit(std::move(fn), std::move(ss), workspace));
    }
    else if (ss->slvrType == SolverType::Implicit) {
      analysis.reset(new FiberRVEAnalysisSImplicit(std::move(fn), std::move(ss), workspace));
    }
    else {
      std::cerr<<"Incorrect Analysis type chosen!"<<std::endl;
      std::abort();
    }
    return analysis;
  }
}
