#ifndef MUMFIM_LINEARTISSUE_H_
#define MUMFIM_LINEARTISSUE_H_
#include <apfFEA.h>
#include <memory>
#include "AnalysisStep.h"
namespace mumfim
{
  class LinearTissueStep : public AnalysisStep
  {
    protected:
    std::map<int, std::unique_ptr<amsi::ElementalSystem> > constitutives;
    public:
    LinearTissueStep(apf::Mesh * mesh,
                 const amsi::ModelDefinition& problem_definition,
                 const amsi::ModelDefinition& solution_strategy,
                 const amsi::ModelDefinition& output,
                 MPI_Comm cm);
    virtual void UpdateDOFs(const double * sol) override;
    virtual void Assemble(amsi::LAS * las) override;
    [[nodiscard]] apf::Field * getUField() final { return apf_primary_field; }
  };
}  // namespace mumfim
#endif
