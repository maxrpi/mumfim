#ifndef MUMFIM_TISSUE_ANALYSIS_H_
#define MUMFIM_TISSUE_ANALYSIS_H_
#include <amsiNonlinearAnalysis.h>
#include "NonlinearTissueStep.h"
#include <model_traits/CategoryNode.h>
namespace mumfim
{

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
    amsi::FEAStep * analysis_step_;
    amsi::LAS * las;
    bool completed;
    // log filenames
    std::string state_fn;
    // logs
    amsi::Log state;
    private:
    static PetscErrorCode CalculateResidual(::SNES s,
                                            Vec solution,
                                            Vec residual,
                                            void * ctx);
    static PetscErrorCode CalculateJacobian(::SNES snes,
                                            Vec displacement,
                                            Mat Amat,
                                            Mat Pmat,
                                            void * ctx);
    static PetscErrorCode CheckConverged(::SNES snes,
                                         PetscInt it,
                                         PetscReal xnorm,
                                         PetscReal gnorm,
                                         PetscReal f,
                                         SNESConvergedReason * reason,
                                         void * ctx);
  };
}  // namespace mumfim
#endif
