#ifndef MUMFIM_TISSUE_ANALYSIS_H_
#define MUMFIM_TISSUE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include "NonlinearTissueStep.h"
#include <model_traits/CategoryNode.h>
namespace mumfim
{
  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param solution_strategy associated solution strategy node
   * @param las the amsi las system used for convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
  template <typename O>
  void buildLASConvergenceOperators(
      const mt::CategoryNode * solution_strategy,
      amsi::MultiIteration * it,
      amsi::LAS * las,
      O out);

  class FEMAnalysis
  {
    public:
    FEMAnalysis(apf::Mesh * mesh,
                   std::unique_ptr<const mt::CategoryNode> cs,
                   MPI_Comm c,
                   const amsi::Analysis & amsi_analysis);
    virtual ~FEMAnalysis();
    virtual void run();
    virtual void checkpoint();
    virtual void finalizeStep();
    virtual void finalizeIteration(int);

    protected:
    // util
    MPI_Comm cm;
    std::unique_ptr<const mt::CategoryNode> analysis_case;
    apf::Mesh* mesh;
    // analysis
    double t;
    double dt;
    int stp;
    int mx_stp;
    int iteration{0};
    AnalysisStep * analysis_step_;
    amsi::LAS * las;
    bool completed;
    // log filenames
    std::string state_fn;
    // logs
    amsi::Log state;
  };
}  // namespace mumfim
#include "TissueAnalysis_impl.h"
#endif
