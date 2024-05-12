#ifndef MUMFIM_NONLINEAR_TISSUE_H_
#define MUMFIM_NONLINEAR_TISSUE_H_
#include "LinearTissueStep.h"
#include "StiffnessVariation.h"
// micro_fo
#include <mumfim/microscale/MultiscaleMicroFOParams.h>
#include <mumfim/microscale/RVE.h>
// amsi
#include <Solvers.h>
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfFunctions.h>
#include <model_traits/CategoryNode.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "amsiFEA.h"

namespace mumfim
{
  class CurrentCoordFunc;
  class NonlinearTissueStep : public amsi::FEAStep
  {
    protected:
    amsi::XpYFunc * xpyfnc;
    std::map<apf::ModelEntity*, std::unique_ptr<amsi::ElementalSystem>> constitutives;
    std::vector<std::unique_ptr<StiffnessVariation>> stf_vrtn_cnst;
    apf::Field * delta_u;
    apf::Field* accepted_displacements;
    apf::Field * current_coords;  // coordinates in current config
    apf::Field * strs;
    apf::Field * strn;
    apf::Field * dfm_grd;
    apf::Field * previous_rve;
    apf::Field * stf_vrtn;
    apf::Field * axl_yngs_mod;
    apf::Field * prev_coords;
    amsi::XpYFunc * prv_crd_fnc;
    double dv_prev;
    int load_step;
    int iteration;
    const mt::CategoryNode & _analysis_case;
    public:
    NonlinearTissueStep(apf::Mesh* mesh, const mt::CategoryNode& analysis_case,
                    MPI_Comm cm = AMSI_COMM_SCALE);
    virtual ~NonlinearTissueStep();
    void computeInitGuess(amsi::LAS* las);
    void step();
    void iter();
    void Assemble(amsi::LAS* las) override;
    void UpdateDOFs(const double* sol) override;
    virtual void preRun() {};
    void AcceptDOFs(){
      apf::copyData(accepted_displacements, getUField());
    }

    int getIteration() { return iteration; }
    apf::Numbering* getNumbering() { return apf_primary_numbering; }
    [[nodiscard]] apf::Field* getUField() final { return apf_primary_field; }
    apf::Mesh* getMesh() { return apf_mesh; }
  };
}
#endif
