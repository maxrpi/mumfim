#ifndef MUMFIM_PROPERTY_EVALUATION_H_
#define MUMFIM_PROPERTY_EVALUATION_H_
#include <amsiNonlinearAnalysis.h>
#include "NonlinearTissueStep.h"
#include <model_traits/CategoryNode.h>
namespace mumfim
{

  class PropertyEvaluation 
  {
    public:
    PropertyEvaluation(apf::Mesh * mesh,
                   std::unique_ptr<const mt::CategoryNode> cs,
                   std::string ktf,
                   MPI_Comm c,
                   const amsi::Analysis & amsi_analysis);
    virtual ~PropertyEvaluation();
    virtual void run();

    protected:
    // util
    MPI_Comm cm;
    std::string kappa_tag_filename;
    std::unique_ptr<const mt::CategoryNode> analysis_case;
    apf::Mesh* mesh;
    amsi::FEAStep * analysis_step_;
  };
}  // namespace mumfim
#endif
