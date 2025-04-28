#include "SinglescaleContinuumAnalysis.h"
#include "StepperFactory.h"
mumfim::SinglescaleContinuumAnalysis::SinglescaleContinuumAnalysis(
    apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
    MPI_Comm c,
    const amsi::Analysis & amsi_analysis,
    std::unordered_map<std::string, std::string> filenames_
  )
    :
    FEMAnalysis(mesh, std::move(cs), c, amsi_analysis)
    , filenames(filenames_)
{
  const auto * solution_strategy =
      mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
  analysis_step_ = createStepper(mesh, *analysis_case, filenames, cm);
}
