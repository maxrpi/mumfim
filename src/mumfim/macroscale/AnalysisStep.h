#ifndef MUMFIM_MACROSCALE_ANALYSISSTEP_H
#define MUMFIM_MACROSCALE_ANALYSISSTEP_H
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>

//#include <iomanip>
//#include <memory>
//#include <unordered_map>

#include "ElementalSystem.h"
//#include "amsiBoundaryConditionQuery.h"
//#include "apfFieldOp.h"
#include "mumfim/analysis/amsiFEA.h"
//#include "mumfim/exceptions.h"

namespace mumfim
{
  class AnalysisStep : public amsi::FEAStep
  {
    protected:

    template <typename T>
    void AssembleIntegratorIntoLAS(amsi::LAS * las,
                                   T integrator,
                                   apf::Field * coordinates = nullptr)
    {
      static_assert(std::is_invocable_r_v<amsi::ElementalSystem *, T,
                                          apf::MeshEntity *, int>);
      if (coordinates == nullptr)
      {
        coordinates = apf_mesh->getCoordinateField();
      }
      apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
      apf::MeshEntity * me = nullptr;
      while ((me = apf_mesh->iterate(it)))
      {
        if (!apf_mesh->isOwned(me))
        {
          continue;
        }
        apf::MeshElement * mlm = apf::createMeshElement(coordinates, me);
        auto * sys = std::invoke(integrator, me, 0);
        apf::Element * elm = apf::createElement(sys->getField(), mlm);
        sys->process(mlm);
        apf::NewArray<apf::Vector3> dofs;
        apf::getVectorNodes(elm, dofs);
        apf::NewArray<int> ids;
        apf::getElementNumbers(apf_primary_numbering, me, ids);
        AssembleDOFs(las, sys->numElementalDOFs(), &ids[0], &dofs[0],
                     &sys->getKe()(0, 0), &sys->getfe()(0),
                     sys->includesBodyForces());
        apf::destroyElement(elm);
        apf::destroyMeshElement(mlm);
      }
      apf_mesh->end(it);
    }

    apf::Mesh * apf_mesh;
    apf::Field * apf_primary_field;
    apf::Field * apf_primary_delta_field;
    apf::Numbering * apf_primary_numbering;
    bool own_mesh;

    public:
    [[nodiscard]] virtual apf::Field * getUField() = 0;
    virtual void preRun() {};
    virtual void step() {};
    virtual void iter() {};
    virtual void AcceptDOFs() {};
    virtual void recoverSecondaryVariables(int) {};
    virtual void computeInitGuess(amsi::LAS * las) {};

    virtual void RenumberDOFs() override;

    AnalysisStep(apf::Mesh * in_mesh,
                 const mt::CategoryNode & analysis_case,
                 std::vector<amsi::DirichletBCEntry> dbc,
                 std::vector<amsi::NeumannBCEntry> nbc,
                 const std::string & analysis_name = "",
                 MPI_Comm cm = AMSI_COMM_SCALE,
                 bool own_mesh = false)
        : FEAStep(analysis_case,
                  std::move(dbc),
                  std::move(nbc),
                  analysis_name,
                  cm)
        , apf_mesh(in_mesh)
        , apf_primary_field(nullptr)
        , apf_primary_delta_field(nullptr)
        , apf_primary_numbering(nullptr)
        , own_mesh(own_mesh)
    {
      analysis_dim = apf_mesh->getDimension();
    }

    AnalysisStep(apf::Mesh * in_mesh,
                 const amsi::ModelDefinition & problem_definition,
                 const amsi::ModelDefinition & solution_strategy,
                 const amsi::ModelDefinition & output,
                 std::vector<amsi::DirichletBCEntry> dbc,
                 std::vector<amsi::NeumannBCEntry> nbc,
                 const std::string & analysis_name = "",
                 MPI_Comm cm = AMSI_COMM_SCALE,
                 bool own_mesh = false)
        : FEAStep(problem_definition,
                  solution_strategy,
                  output,
                  std::move(dbc),
                  std::move(nbc),
                  analysis_name,
                  cm)
        , apf_mesh(in_mesh)
        , apf_primary_field(nullptr)
        , apf_primary_delta_field(nullptr)
        , apf_primary_numbering(nullptr)
        , own_mesh(own_mesh)
    {
      analysis_dim = apf_mesh->getDimension();
    }

    ~AnalysisStep()
    {
      if (own_mesh)
      {
        apf::destroyMesh(apf_mesh);
      }
    }

    apf::Mesh * getMesh() { return apf_mesh; }

    virtual void UpdateDOFs(const double *) override;
    virtual void ApplyBC_Dirichlet() override;
    virtual void ApplyBC_Neumann(amsi::LAS * las) override;
    virtual void Adapt() override;
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_MULTISCALETISSUE_CC_TISSUEBASE_H
