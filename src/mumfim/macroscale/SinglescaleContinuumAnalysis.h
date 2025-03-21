#ifndef MUMFIM_SINGLESCALECONTINUUMANALYSIS_H
#define MUMFIM_SINGLESCALECONTINUUMANALYSIS_H
#include "FEMAnalysis.h"

namespace mumfim
{
  class SinglescaleContinuumAnalysis : public FEMAnalysis
  {
    public:
    SinglescaleContinuumAnalysis(apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
        MPI_Comm c,
    const amsi::Analysis & amsi_analysis);
  };
}  // namespace mumfim
#endif
