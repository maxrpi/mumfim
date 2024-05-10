#include "MultiscaleTissueAnalysis.h"
#include <Solvers.h>
#include <amsiCasters.h>
#include <amsiMultiscale.h>
#include <apfFunctions.h>
#include <apfNumbering.h>
#include <apfWrapper.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include "ModelTraits.h"
#include "MultiscaleTissueStep.h"
#include "mumfim/exceptions.h"
namespace mumfim
{
  void MultiscaleTissueIteration::iterate()
  {
    std::cerr<<"We shouldn't be calling the multiscale tissue iteration if code is working as expected\n";
    std::abort();
  }
  MultiscaleTissueAnalysis::MultiscaleTissueAnalysis(
      apf::Mesh * mesh,
      std::unique_ptr<mt::CategoryNode> analysis_case,
      MPI_Comm cm,
      const amsi::Analysis & amsi_analysis,
      const amsi::Multiscale & amsi_multiscale)
      : FEMAnalysis(mesh, std::move(analysis_case), cm, amsi_analysis)
      , cplng(getRelationID(amsi_multiscale.getMultiscaleManager(),
                            amsi_multiscale.getScaleManager(),
                            "macro",
                            "micro_fo"))
      , multiscale_(amsi_multiscale)
  {
    const auto * solution_strategy =
        mt::GetPrimaryCategoryByType((this->analysis_case.get()), "solution strategy");
    if(solution_strategy == nullptr) {
      throw mumfim_error("solution strategy is required for multiscale analysis");
    }
    analysis_step_ = new MultiscaleTissueStep(mesh, *(this->analysis_case), cm, multiscale_);
    static_cast<MultiscaleTissueStep *>(analysis_step_)->initMicro();
  }
  void MultiscaleTissueAnalysis::finalizeStep()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(cplng, &completed);
  }

  void MultiscaleTissueAnalysis::finalizeIteration(int lastIteration) {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(cplng, &lastIteration);
  }
}  // namespace mumfim
