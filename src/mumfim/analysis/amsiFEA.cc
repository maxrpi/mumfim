#include "amsiFEA.h"

#include <apfNumbering.h>
#include <apfShape.h>
#include <maShape.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <vector>

#include "amsiBoundaryConditions.h"
#include "apfFunctions.h"

namespace amsi
{
  void FEAStep::RenumberDOFs()
  {
    if (!numbered)
    {
      local_dof_count = apf::naiveOrder(apf_primary_numbering);
      int analysis_size, analysis_rank;
      MPI_Comm_rank(analysis_comm, &analysis_rank);
      MPI_Comm_size(analysis_comm, &analysis_size);
      // add constraint DOFs to the last rank
      if (analysis_size - 1 == analysis_rank)
        local_dof_count += constraint_dofs;
      // FIXME use nonblocking calls here.
      MPI_Exscan(&local_dof_count, &first_local_dof, 1, MPI_INTEGER, MPI_SUM,
                 analysis_comm);
      // need to set the global dof counts for the fea class
      // which is used externally.
      MPI_Allreduce(&local_dof_count, &global_dof_count, 1, MPI_INTEGER,
                    MPI_SUM, analysis_comm);
      apf::SetNumberingOffset(apf_primary_numbering, first_local_dof);
      apf::synchronize(apf_primary_numbering);
      numbered = true;
    }
  }  // use solution vector to update displacement dofs associated with

  // locally-owned nodes
  void FEAStep::UpdateDOFs(const double * solution)
  {
    int num_components = apf::countComponents(apf_primary_field);
    apf::MeshEntity * mesh_ent = NULL;
    for (int ii = 0; ii < analysis_dim; ii++)
    {
      for (apf::MeshIterator * it = apf_mesh->begin(ii);
           (mesh_ent = apf_mesh->iterate(it));)
      {
        if (apf_mesh->isOwned(mesh_ent))
        {
          apf::FieldShape * fs = apf::getShape(apf_primary_field);
          int num_nodes = fs->countNodesOn(apf_mesh->getType(mesh_ent));
          for (int jj = 0; jj < num_nodes; jj++)
          {
            apf::Vector3 disp;
            apf::getVector(apf_primary_field, mesh_ent, jj, disp);
            for (int kk = 0; kk < num_components; kk++)
            {
              if (!apf::isFixed(apf_primary_numbering, mesh_ent, jj, kk))
              {
                int global_number =
                    getNumber(apf_primary_numbering, mesh_ent, jj, kk);
                disp[kk] += solution[global_number - first_local_dof];
              }
            }
            apf::setVector(apf_primary_field, mesh_ent, jj, disp);
          }
        }
      }
    }
    apf::synchronize(apf_primary_field);
  }

  void FEAStep::ApplyBC_Dirichlet()
  {
    fixed_dofs =
        applyDirichletBCs(apf_primary_numbering, apf_primary_delta_field,
                          *problem_definition.associated, dirichlet_bcs, T);
    fixed_dofs = comm_sum(fixed_dofs, analysis_comm);
    AMSI_V1(std::cout << "There are " << fixed_dofs
                      << " dofs fixed by essential boundary conditions."
                      << std::endl;)
  }

  void FEAStep::ApplyBC_Neumann(amsi::LAS * las)
  {
    applyNeumannBCs(las, apf_primary_numbering, *problem_definition.associated,
                    neumann_bcs, T);
  }

  void FEAStep::Adapt()
  {
    throw mumfim::mumfim_error("Adapt Loop not implemented");
  }

  FEAStep::FEAStep(apf::Mesh * mesh,
                   const ModelDefinition & pd,
                   const ModelDefinition & ss,
                   const ModelDefinition & out,
                   std::vector<DirichletBCEntry> dbc,
                   std::vector<NeumannBCEntry> nbc,
                   const std::string & analysis_name,
                   MPI_Comm cm,
                   bool own_mesh)
      : constraint_dofs(0)
      , local_dof_count(0)
      , first_local_dof(0)
      , global_dof_count(0)
      , fixed_dofs(0)
      , numbered(false)
      , T(0.0)
      , analysis_comm(cm)
      , problem_definition(pd)
      , output(out)
      , solution_strategy(ss)
      , dirichlet_bcs(std::move(dbc))
      , neumann_bcs(std::move(nbc))
      , apf_mesh(mesh)
      , apf_primary_field(nullptr)
      , apf_primary_delta_field(nullptr)
      , apf_primary_numbering(nullptr)
      , own_mesh(own_mesh)
  {
    analysis_dim = apf_mesh->getDimension();
  }

  static ModelDefinition GetModelDefinitionComponent(
      const mt::CategoryNode & analysis_case,
      const std::string & analysis_name, const std::string& component_name)
  {
    auto* component = mt::GetCategoryByType(&analysis_case, component_name);
    auto* temp = mt::GetCategoryByType(component, analysis_name);
    component = (temp == nullptr) ? component : temp;
    if(component == nullptr) {
      throw mumfim::mumfim_error(component_name + " not defined");
    }
    return {
        .associated = mt::AssociateModel(component),
        .unassociated = std::make_shared<mt::CategoryNode>(*component)};
  }

  FEAStep::FEAStep(apf::Mesh * mesh,
                   const mt::CategoryNode & analysis_case,
                   std::vector<DirichletBCEntry> dbc,
                   std::vector<NeumannBCEntry> nbc,
                   const std::string & analysis_name,
                   MPI_Comm cm,
                   bool own_mesh)
      : FEAStep(
            mesh,
            GetModelDefinitionComponent(analysis_case,
                                        analysis_name,
                                        "problem definition"),
            GetModelDefinitionComponent(analysis_case,
                                        analysis_name,
                                        "solution strategy"),
            GetModelDefinitionComponent(analysis_case, analysis_name, "output"),
            std::move(dbc),
            std::move(nbc),
            analysis_name,
            cm,
            own_mesh)
  {
  }

  FEAStep::~FEAStep()
  {
    if (own_mesh)
    {
      apf::destroyMesh(apf_mesh);
    }
  }

  void FEAStep::setSimulationTime(double t)
  {
    T = t;
    if (!PCU_Comm_Self())
      std::cout << "Simulation time updated: " << T << std::endl;
  }

  void FEAStep::GetDOFInfo(int & global, int & local, int & offset)
  {
    global = global_dof_count;
    local = local_dof_count;
    offset = first_local_dof;
  }

  template <typename NODE_TYPE>
  static void AssembleDOFs(LAS* las, int num_elemental_dofs, int* dof_numbers,
                             const NODE_TYPE* node_values, double* Ke, double* fe,
                             bool includes_body_forces, int analysis_dim)
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
  void AssembleDOFs<double>(LAS* las, int num_elemental_dofs, int* dof_numbers,
                                            const double* node_values, double* Ke, double* fe,
                                            bool includes_body_forces, int analysis_dim)
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


  void FEAStep::AssembleIntegratorIntoLAS(LAS * las,
                                 apf::Field * coordinates)

  {
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
      auto* sys = getIntegrator(me, 0);
      apf::Element * elm = apf::createElement(sys->getField(), mlm);
      sys->process(mlm);
      apf::NewArray<apf::Vector3> dofs;
      apf::getVectorNodes(elm, dofs);
      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering, me, ids);
      AssembleDOFs(las, sys->numElementalDOFs(), &ids[0], &dofs[0],
                   &sys->getKe()(0, 0), &sys->getfe()(0),
                   sys->includesBodyForces(), analysis_dim);
      apf::destroyElement(elm);
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
  }

}  // namespace amsi
