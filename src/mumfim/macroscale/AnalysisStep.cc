#include "AnalysisStep.h"

#include <apfNumbering.h>
#include <apfShape.h>
#include <maShape.h>

#include <cassert>

#include "amsiBoundaryConditions.h"
#include "apfFunctions.h"

namespace mumfim {
void AnalysisStep::RenumberDOFs()
{
  if (!numbered) {
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
}
  // use solution vector to update displacement dofs associated with
  // locally-owned nodes. This assumes 'solution' is an increment,
  // i.e., a delta, on the apf_primary_field.
  void AnalysisStep::UpdateDOFs(const double* solution)
  {
    if(apf::getValueType(apf_primary_field) == apf::VECTOR){
      int num_components = apf::countComponents(apf_primary_field);
      apf::MeshEntity* mesh_ent = NULL;
      for (int ii = 0; ii < analysis_dim; ii++) {
        for (apf::MeshIterator* it = apf_mesh->begin(ii);
            (mesh_ent = apf_mesh->iterate(it));) {
          if (apf_mesh->isOwned(mesh_ent)) {
            apf::FieldShape* fs = apf::getShape(apf_primary_field);
            int num_nodes = fs->countNodesOn(apf_mesh->getType(mesh_ent));
            for (int jj = 0; jj < num_nodes; jj++) {
              apf::Vector3 disp;
              apf::getVector(apf_primary_field, mesh_ent, jj, disp);
              for (int kk = 0; kk < num_components; kk++) {
                if (!apf::isFixed(apf_primary_numbering, mesh_ent, jj, kk)) {
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
    } else {
      int num_components = apf::countComponents(apf_primary_field);
      apf::MeshEntity* mesh_ent = NULL;
      for (int ii = 0; ii < analysis_dim; ii++) {
        for (apf::MeshIterator* it = apf_mesh->begin(ii);
            (mesh_ent = apf_mesh->iterate(it));) {
          if (apf_mesh->isOwned(mesh_ent)) {
            apf::FieldShape* fs = apf::getShape(apf_primary_field);
            int num_nodes = fs->countNodesOn(apf_mesh->getType(mesh_ent));
            for (int jj = 0; jj < num_nodes; jj++) {
              double disp = apf::getScalar(apf_primary_field, mesh_ent, jj);
              if (!apf::isFixed(apf_primary_numbering, mesh_ent, jj, 0)) {
                int global_number =
                    getNumber(apf_primary_numbering, mesh_ent, jj, 0);
                disp += solution[global_number - first_local_dof];
              }
              apf::setScalar(apf_primary_field, mesh_ent, jj, disp);
            }
          }
        }
      }
    }
    apf::synchronize(apf_primary_field);
  }  // FIXME Apply boundary conditions without simmetrix
  void AnalysisStep::ApplyBC_Dirichlet()
  {
    fixed_dofs =
        applyDirichletBCs(apf_primary_numbering, apf_primary_delta_field,
                          *problem_definition.associated, dirichlet_bcs, T);
    fixed_dofs = amsi::comm_sum(fixed_dofs, analysis_comm);
    AMSI_V1(std::cout << "There are " << fixed_dofs
                      << " dofs fixed by essential boundary conditions."
                      << std::endl;)
  }  // FIXME apply boundary conditions without simmetrix
  void AnalysisStep::ApplyBC_Neumann(amsi::LAS* las)
  {
    applyNeumannBCs(las, apf_primary_numbering, *problem_definition.associated,
                    neumann_bcs, T);
  }
  void AnalysisStep::Adapt()
  {
    throw mumfim_error("Adapt Loop not implemented");
  }
}
