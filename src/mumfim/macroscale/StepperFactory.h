#ifndef MUMFIM_STEPPER_FACTORY_H_
#define MUMFIM_STEPPER_FACTORY_H_
#include <amsiNonlinearAnalysis.h>
#include <amsiMultiscale.h>
#include "amsiFEA.h"
#include <unordered_map>

namespace mumfim
{
   amsi::FEAStep* createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      std::unordered_map<std::string, std::string> filenames,
      MPI_Comm comm
    );

    amsi::FEAStep * createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      const amsi::Multiscale & amsi_multiscale,
      MPI_Comm comm
    );

}  // namespace mumfim;
#endif