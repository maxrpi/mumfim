#ifndef MUMFIM_TISSUEMULTISCALEANALYSIS_H_
#define MUMFIM_TISSUEMULTISCALEANALYSIS_H_
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiPETScLAS.h>
#include <apf.h>
#include <model_traits/CategoryNode.h>
#include <stdexcept>
#include <string>
#include "LinearTissueStep.h"
#include "MultiscaleTissueStep.h"
#include "FEMAnalysis.h"
#include "VolumeConvergence.h"
namespace mumfim
{
  class MultiscaleTissueIteration : public amsi::Iteration
  {
    protected:
    MultiscaleTissueStep * tssu;
    amsi::LAS * las;
    Iteration * fem_iter;

    public:
    MultiscaleTissueIteration(MultiscaleTissueStep * a, amsi::LAS * l)
        : tssu(a), las(l), fem_iter(amsi::buildLinearFEMIteration(a, l))
    {
    }
    ~MultiscaleTissueIteration() { delete fem_iter; }
    void iterate();
  };
  class MultiscaleTissueAnalysis : public FEMAnalysis
  {
    public:
    MultiscaleTissueAnalysis(apf::Mesh * mesh,
                             std::unique_ptr<mt::CategoryNode> analysis_case,
                             MPI_Comm cm,
                             const amsi::Analysis & amsi_analysis,
                             const amsi::Multiscale & amsi_multiscale);
    virtual void finalizeStep() final;
    virtual void finalizeIteration(int) final;


    private:
    size_t cplng;
    const amsi::Multiscale & multiscale_;
  };
}  // namespace mumfim
#endif
