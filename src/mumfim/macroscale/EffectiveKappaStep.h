#ifndef MUMFIM_LINEAR_HEAT_CONDUCTION_STEP_H_
#define MUMFIM_LINEAR_HEAT_CONDUCTION_STEP_H_
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
  class EffectiveKappaStep : public amsi::FEAStep
  {
    protected:
    std::map<apf::ModelEntity*, std::unique_ptr<amsi::ElementalSystem>> constitutives;
    apf::Field * temperature;
    apf::Field * flux;
    apf::Field * kappa;
    apf::Field * sources;
    std::map<std::pair<int, int>, double *> k_map;
    int iteration;
    std::vector <int> interior, exterior; 
    double RVEvolume;

    [[nodiscard]] amsi::ElementalSystem * getIntegrator(
        apf::MeshEntity * ent,
        int integration_point) override;

    public:
    EffectiveKappaStep(apf::Mesh* mesh, const mt::CategoryNode& analysis_case,
                    MPI_Comm comm_ = AMSI_COMM_SCALE);
    virtual ~EffectiveKappaStep();
    void computeInitGuess(amsi::LAS* las);
    void Assemble(amsi::LAS* las) override;
    apf::Matrix3x3 computeKM(const mt::CategoryNode& analysis_case, MPI_Comm comm_ = AMSI_COMM_SCALE);
    apf::Field* getUField() { return apf_primary_field; }

  };


}
#endif