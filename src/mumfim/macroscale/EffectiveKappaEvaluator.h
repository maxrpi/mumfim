#ifndef MUMFIM_EFFECTIVE_KAPPA_EVALUATOR_H_
#define MUMFIM_EFFECTIVE_KAPPA_EVALUATOR_H_
// amsi
#include <Solvers.h>
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfFunctions.h>
#include <model_traits/CategoryNode.h>
#include <petscmat.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "amsiFEA.h"

namespace mumfim
{
  class EffectiveKappaEvaluator : public amsi::FEAStep
  {
    protected:
    std::map<apf::ModelEntity*, std::unique_ptr<amsi::ElementalSystem>> constitutives;
    apf::Field * kappa;
    apf::Field * coordinates;
    int iteration;
    apf::Matrix3x3 Km;
    // Mesh vertex number to partitioned vertex number
    std::vector<bool> onExterior;
    std::vector<int> vert2subvert;
    apf::Vector3 centroid; //centroid of boundary nodes
    int n_int = 0;
    int n_ext = 0;
    Mat Kii, Kie, Kei, Kee;

    [[nodiscard]] amsi::ElementalSystem * getIntegrator(
        apf::MeshEntity * ent,
        int integration_point) override;

    void collectInteriorAndExteriorVertices(void);
    void AssembleDOFs(const std::vector<int> & dof_numbers,
                         apf::DynamicMatrix & Ke);
    void locateVertices(void);
    void AssembleIntegratorIntoMat(void);

    public:
    EffectiveKappaEvaluator(apf::Mesh* mesh, const mt::CategoryNode& analysis_case,
                    MPI_Comm comm_ = AMSI_COMM_SCALE);
    virtual ~EffectiveKappaEvaluator();
    void Assemble(amsi::LAS*);
    void Solve();
    apf::Field* getUField() { return apf_primary_field; }

  };


}
#endif