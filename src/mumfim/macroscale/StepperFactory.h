#ifndef MUMFIM_STEPPER_FACTORY_H_
#define MUMFIM_STEPPER_FACTORY_H_
#include <amsiNonlinearAnalysis.h>
#include <amsiMultiscale.h>
#include "amsiFEA.h"
namespace mumfim
{
   amsi::FEAStep* createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      MPI_Comm comm,
      std::string ktf = ""
    );

    amsi::FEAStep * createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      const amsi::Multiscale & amsi_multiscale,
      MPI_Comm comm
    );

}  // namespace mumfim;
#endif