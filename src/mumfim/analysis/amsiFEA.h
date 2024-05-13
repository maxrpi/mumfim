#ifndef AMSI_FEA_H_
#define AMSI_FEA_H_
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <mpi.h>

#include <cstring>  // memset
#include <iostream>
#include <list>
#include <map>
#include <vector>

#include "ElementalSystem.h"
#include "amsiBoundaryConditions.h"
#include "amsiLAS.h"
#include "amsiMPI.h"
#include "model_traits/AssociatedModelTraits.h"
#include "mumfim/analysis/amsiFEA.h"
#include "mumfim/exceptions.h"
namespace amsi {
  struct ModelDefinition {
    std::shared_ptr<const mt::AssociatedModelTraits<mt::DimIdGeometry>> associated;
    std::shared_ptr<const mt::CategoryNode> unassociated;
  };
  typedef int ldof_type;
  typedef long gdof_type;
  // TODO: replace with functor subclassed on field library... or something
  bool isFixed(int);
  template <typename T>
  void Model_PrintInfo(T model, std::ostream& out);
  template <typename T>
  void Mesh_PrintInfo(T mesh, std::ostream& out);
  class FEAStep
  {
    protected:
    // additional global dofs from constraints
    int constraint_dofs;
    // the number of dofs on this processor
    int local_dof_count;
    // the first locally-owned dof (using the global numbering)
    int first_local_dof;
    // the total number of dofs in the analysis
    int global_dof_count;
    // the number of fixed locally-fixed dofs
    int fixed_dofs;
    // whether the dofs have been numbered or not
    bool numbered;
    // the current simulation time (for dynamics) or psuedo-time (for statics)
    // might want to pull this out entirely somehow... just pass in to pertinent
    // member functions especially as it's really only used for nonlinear FEAStep,
    // for linear FEAStep it is just set to 1
    double T;
    // name of this analysis
    std::string name;
    // dimensionality of the analysis
    int analysis_dim;
    // mpi communicator containing only those processes undertaking this
    // analysis
    MPI_Comm analysis_comm;
    ModelDefinition problem_definition;
    ModelDefinition output;
    ModelDefinition solution_strategy;
    std::vector<DirichletBCEntry> dirichlet_bcs;
    std::vector<NeumannBCEntry> neumann_bcs;
    /**
    * this function is a hook for subclasses to specify what elemental system
    * should be used within each mesh entity
    * \param ent current mesh entity
    * \param integration_point current integration point
    */
    [[nodiscard]] virtual amsi::ElementalSystem * getIntegrator(
        apf::MeshEntity * ent,
        int integration_point) = 0;

    void AssembleIntegratorIntoLAS(LAS * las,
                                   apf::Field * coordinates = nullptr);

    apf::Mesh * apf_mesh;
    apf::Field * apf_primary_field;
    apf::Field * apf_primary_delta_field;
    apf::Numbering * apf_primary_numbering;
    bool own_mesh;

    public:
    FEAStep(apf::Mesh * mesh,
            const mt::CategoryNode & analysis_case,
            std::vector<DirichletBCEntry> dbc,
            std::vector<NeumannBCEntry> nbc,
            const std::string & analysis_name = "",
            MPI_Comm cm = AMSI_COMM_SCALE,
            bool own_mesh = false);

    FEAStep(apf::Mesh * mesh,
            const ModelDefinition & problem_definition,
            const ModelDefinition & solution_strategy,
            const ModelDefinition & output,
            std::vector<DirichletBCEntry> dbc,
            std::vector<NeumannBCEntry> nbc,
            const std::string & analysis_name = "",
            MPI_Comm cm = AMSI_COMM_SCALE,
            bool own_mesh = false);
    virtual ~FEAStep();
    virtual void Assemble(LAS*) = 0;
    void setSimulationTime(double t);
    virtual void GetDOFInfo(int& global, int& local, int& offset);

    virtual void RenumberDOFs();
    virtual void UpdateDOFs(const double *);
    [[nodiscard]] virtual apf::Field * getUField() = 0;

    virtual void preRun() {}

    virtual void step() {}

    virtual void iter() {}

    virtual void AcceptDOFs() {}

    virtual void recoverSecondaryVariables(int) {}

    virtual void computeInitGuess(LAS * las) {}

    virtual void ApplyBC_Dirichlet();
    virtual void ApplyBC_Neumann(LAS * las);
    virtual void Adapt();

    apf::Mesh * getMesh() { return apf_mesh; }
  };

}  // namespace amsi
#endif
