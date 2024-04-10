#ifndef MUMFIM_SINGLESCALETISSUEANALYSIS_H
#define MUMFIM_SINGLESCALETISSUEANALYSIS_H
#include "FEMAnalysis.h"

namespace mumfim
{
  class SinglescaleTissueAnalysis : public FEMAnalysis
  {
    public:
    SinglescaleTissueAnalysis(apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
        MPI_Comm c,
    const amsi::Analysis & amsi_analysis);
  };
}  // namespace mumfim
#endif
