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

    template <typename T>
    void AssembleIntegratorIntoLAS(LAS * las,
                                   T integrator,
                                   apf::Field * coordinates = nullptr)
    {
      static_assert(std::is_invocable_r_v<ElementalSystem *, T,
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
    template <typename NODE_TYPE>
    void AssembleDOFs(LAS* las, int num_elemental_dofs, int* dof_numbers,
                      const NODE_TYPE* node_values, double* Ke, double* fe,
                      bool includes_body_forces) const;
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

  template <typename NODE_TYPE>
  void FEAStep::AssembleDOFs(LAS* las, int num_elemental_dofs, int* dof_numbers,
                         const NODE_TYPE* node_values, double* Ke, double* fe,
                         bool includes_body_forces) const
  {
    if (Ke != NULL) {
      if (!includes_body_forces) {
        double* bf = new double[num_elemental_dofs]();
        for (int ii = 0; ii < num_elemental_dofs; ii++) {
          const int& global_i = dof_numbers[ii];
          // this is the isFixed function from apfFunctions
          // which is different from the sim query is fixed function!
          if (!isFixed(global_i)) {
            for (int jj = 0; jj < num_elemental_dofs; jj++) {
              const double& val = Ke[ii * num_elemental_dofs + jj];
              double j_val = node_values[jj / analysis_dim][jj % analysis_dim];
              if (j_val != 0.0) bf[ii] += -val * j_val;
            }
          }
        }
        las->AddToVector(num_elemental_dofs, dof_numbers, &bf[0]);
        delete[] bf;
      }
      las->AddToMatrix(num_elemental_dofs, dof_numbers, num_elemental_dofs,
                       dof_numbers, &Ke[0]);
    }
    /// Modification of fe to correctly account for nonzero dirichlet boundary
    /// conditions
    double* dirichletValue = new double[num_elemental_dofs]();
    for (int ii = 0; ii < num_elemental_dofs; ii++) {
      if (dof_numbers[ii] < 0)
        dirichletValue[ii] = node_values[ii / analysis_dim][ii % analysis_dim];
    }
    double* dfe = new double[num_elemental_dofs]();
    for (int ii = 0; ii < num_elemental_dofs; ii++)
      for (int jj = 0; jj < num_elemental_dofs; jj++)
        dfe[ii] =
            dfe[ii] + Ke[ii * num_elemental_dofs + jj] * dirichletValue[ii];
    for (int ii = 0; ii < num_elemental_dofs; ii++)
      fe[ii] = fe[ii] + dfe[ii];
    delete[] dirichletValue;
    delete[] dfe;
    las->AddToVector(num_elemental_dofs, dof_numbers, &fe[0]);
  }

  template <>
  inline void FEAStep::AssembleDOFs<double>(LAS* las, int num_elemental_dofs, int* dof_numbers,
                         const double* node_values, double* Ke, double* fe,
                         bool includes_body_forces) const
  {
    if (Ke == NULL) {
      throw mumfim::mumfim_error("Ke is null going into AssembleDOFs");
    }
    if (!includes_body_forces) {
      throw mumfim::mumfim_error("Scalar does not handle case with !includes_body_force");
    }

    las->AddToMatrix(num_elemental_dofs, dof_numbers, num_elemental_dofs,
                     dof_numbers, &Ke[0]);

    /// Modification of fe to correctly account for nonzero dirichlet boundary
    /// conditions
    std::vector<double> dirichletValue(num_elemental_dofs, 0.0);
    for (int ii = 0; ii < num_elemental_dofs; ii++) {
      if (dof_numbers[ii] < 0) dirichletValue[ii] = node_values[ii];
    }
    std::vector<double> dfe(num_elemental_dofs, 0.0);
    for (int ii = 0; ii < num_elemental_dofs; ii++)
      for (int jj = 0; jj < num_elemental_dofs; jj++)
        dfe[ii] =
            dfe[ii] + Ke[ii * num_elemental_dofs + jj] * dirichletValue[ii];
    for (int ii = 0; ii < num_elemental_dofs; ii++)
      fe[ii] = fe[ii] + dfe[ii];
    las->AddToVector(num_elemental_dofs, dof_numbers, &fe[0]);
  }
}  // namespace amsi
#endif
