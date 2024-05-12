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
  class LinearHeatConductionStep : public amsi::FEAStep
  {
    protected:
    std::map<apf::ModelEntity*, std::unique_ptr<amsi::ElementalSystem>> constitutives;
    apf::Field * temperature;
    apf::Field * flux;
    apf::Field * kappa;
    apf::Field * sources;
    int iteration;
    public:
    LinearHeatConductionStep(apf::Mesh* mesh, const mt::CategoryNode& analysis_case,
                    MPI_Comm comm_ = AMSI_COMM_SCALE);
    virtual ~LinearHeatConductionStep();
    void computeInitGuess(amsi::LAS* las);
    //void step();
    //void iter();
    void Assemble(amsi::LAS* las) override;
    //void UpdateDOFs(const double* sol) override;
    //virtual void recoverSecondaryVariables(int);
    //virtual void preRun() {};
    //void AcceptDOFs(){ apf::copyData(temperature, getUField()); }
    //int getIteration() { return iteration; }
    //apf::Numbering* getNumbering() { return apf_primary_numbering; }
    apf::Field* getUField() { return apf_primary_field; }
    //apf::Mesh* getMesh() { return apf_mesh; }
    //void getFluxOn(apf::ModelEntity* ent, double* frc);
  };


}
#endif