#include "amsiFEA.h"

#include <apfFieldData.h>
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
          apf::FieldShape * shape_function = apf::getShape(apf_primary_field);
          int num_nodes =
              shape_function->countNodesOn(apf_mesh->getType(mesh_ent));
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
      const std::string & analysis_name,
      const std::string & component_name)
  {
    auto * component = mt::GetCategoryByType(&analysis_case, component_name);
    auto * temp = mt::GetCategoryByType(component, analysis_name);
    component = (temp == nullptr) ? component : temp;
    if (component == nullptr)
    {
      throw mumfim::mumfim_error(component_name + " not defined");
    }
    return {.associated = mt::AssociateModel(component),
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

  static void AssembleDOFs(LAS * las,
                           const std::vector<int> & dof_numbers,
                           const std::vector<double> & dof_values,
                           apf::DynamicMatrix & Ke,
                           apf::DynamicVector & fe,
                           bool includes_body_forces,
                           int ncomponents)
  {
    if (!includes_body_forces)
    {
      throw mumfim::mumfim_error(
          "Integrator that does not bring int it's own body forces is no "
          "longer supported");
    }
    const auto num_elemental_dofs = fe.size();
    las->AddToMatrix(num_elemental_dofs, dof_numbers.data(), num_elemental_dofs,
                     dof_numbers.data(), &Ke(0, 0));
    for (int ii = 0; ii < num_elemental_dofs; ii++)
    {
      // if the dof_number is negative than that particular dof is fixed
      // so set the value to zero if it's not fixed and to the real value if
      // it is fixed.
      double dirichlet_value = dof_numbers[ii] < 0 ? dof_values[ii] : 0.0;
      for (int jj = 0; jj < num_elemental_dofs; jj++)
      {
        fe[ii] += Ke(ii, jj) * dirichlet_value;
      }
    }
    las->AddToVector(num_elemental_dofs, dof_numbers.data(), &fe[0]);
  }

  void FEAStep::AssembleIntegratorIntoLAS(LAS * las, apf::Field * coordinates)

  {
    if (coordinates == nullptr)
    {
      coordinates = apf_mesh->getCoordinateField();
    }
    apf::MeshIterator * mesh_region_iter = apf_mesh->begin(analysis_dim);
    while (apf::MeshEntity * mesh_entity = apf_mesh->iterate(mesh_region_iter))
    {
      if (!apf_mesh->isOwned(mesh_entity))
      {
        continue;
      }
      apf::MeshElement * mlm = apf::createMeshElement(coordinates, mesh_entity);
      auto * sys = getIntegrator(mesh_entity, 0);
      sys->process(mlm);

      AssembleDOFs(las, sys->getFieldNumbers(), sys->getFieldValues(),
                   sys->getKe(), sys->getfe(), sys->includesBodyForces(),
                   apf::countComponents(sys->getField()));
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(mesh_region_iter);
  }

}  // namespace amsi
