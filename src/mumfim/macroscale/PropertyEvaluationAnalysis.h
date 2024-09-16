#ifndef MUMFIM_PROPERTYEVALUATIONANALYSIS_H
#define MUMFIM_PROPERTYEVALUATIONANALYSIS_H
#include "FEMAnalysis.h"

namespace mumfim
{
  class PropertyEvaluationAnalysis : public FEMAnalysis
  {
    public:
    PropertyEvaluationAnalysis(apf::Mesh * mesh,
    std::unique_ptr<const mt::CategoryNode> cs,
        MPI_Comm c,
    const amsi::Analysis & amsi_analysis);
  };
}  // namespace mumfim
#endif
