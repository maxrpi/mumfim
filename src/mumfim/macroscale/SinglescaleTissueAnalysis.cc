#include "SinglescaleTissueAnalysis.h"
#include "StepperFactory.h"
mumfim::SinglescaleTissueAnalysis::SinglescaleTissueAnalysis(
    apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
    MPI_Comm c,
    const amsi::Analysis & amsi_analysis)
    : FEMAnalysis(mesh, std::move(cs), c, amsi_analysis)
{
  const auto * solution_strategy =
      mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
  analysis_step_ = createStepper(mesh, *analysis_case, cm);
}
