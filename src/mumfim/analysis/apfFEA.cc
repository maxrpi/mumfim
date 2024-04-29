#include "apfFEA.h"
#include <apfNumbering.h>
#include <apfShape.h>
#include <maShape.h>
#include <cassert>
#include "amsiBoundaryConditions.h"
#include "apfFunctions.h"
namespace amsi {
  // assumes only linear tets
  apf::Field* analyzeMeshQuality(apf::Mesh* mesh, apf::Field* disp_field)
  {
    int analysis_dim = mesh->getDimension();
    apf::Field* f = createStepField(mesh, "quality", apf::SCALAR);
    apf::MeshEntity* me = NULL;
    for (apf::MeshIterator* it = mesh->begin(analysis_dim);
         (me = mesh->iterate(it));) {
      apf::Vector3 verts[4];
      apf::Adjacent v;
      mesh->getAdjacent(me, 0, v);
      for (int ii = 0; ii < 4; ii++) {
        apf::Vector3 disp;
        apf::getVector(disp_field, v[ii], 0, disp);
        mesh->getPoint(v[ii], 0, verts[ii]);
        verts[ii] += disp;
      }
      apf::setScalar(f, me, 0, ma::measureLinearTetQuality(verts));
    }
    return f;
  }
  void apfFEA::RenumberDOFs()
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
  void apfFEA::Assemble(LAS* las)
  {
    assert(elemental_system);
    apf::MeshEntity* me = NULL;
    apf::MeshIterator* it = apf_mesh->begin(analysis_dim);
    while ((me = apf_mesh->iterate(it))) {
      apf::MeshElement* melm = apf::createMeshElement(apf_mesh, me);
      apf::Element* elm =
          apf::createElement(elemental_system->getField(), melm);
      elemental_system->process(melm);

      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering, me, ids);

      if( apf::countComponents(elemental_system->getField()) == 1)
      {
        apf::NewArray<double> dofs;
        apf::getScalarNodes(elm, dofs);
        AssembleDOFs(las, elemental_system->numElementalDOFs(), &ids[0], &dofs[0],
                   &elemental_system->getKe()(0, 0),
                   &elemental_system->getfe()(0),
                   elemental_system->includesBodyForces());
      } else 
      {
        apf::NewArray<apf::Vector3> dofs;
        apf::getVectorNodes(elm, dofs);
        AssembleDOFs(las, elemental_system->numElementalDOFs(), &ids[0], &dofs[0],
                   &elemental_system->getKe()(0, 0),
                   &elemental_system->getfe()(0),
                   elemental_system->includesBodyForces());
      }
      apf::destroyElement(elm);
      apf::destroyMeshElement(melm);
    }
    apf_mesh->end(it);
  }
  // use solution vector to update displacement dofs associated with
  // locally-owned nodes
  void apfFEA::UpdateDOFs(const double* solution)
  {
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
    apf::synchronize(apf_primary_field);
  }
  // FIXME Apply boundary conditions without simmetrix
  void apfFEA::ApplyBC_Dirichlet()
  {
    fixed_dofs =
        applyDirichletBCs(apf_primary_numbering, apf_primary_delta_field,
                          *problem_definition.associated, dirichlet_bcs, T);
    fixed_dofs = comm_sum(fixed_dofs, analysis_comm);
    AMSI_V1(std::cout << "There are " << fixed_dofs
                      << " dofs fixed by essential boundary conditions."
                      << std::endl;)
  }
  // FIXME apply boundary conditions without simmetrix
  void apfFEA::ApplyBC_Neumann(LAS* las)
  {
    applyNeumannBCs(las, apf_primary_numbering, *problem_definition.associated,
                    neumann_bcs, T);
  }
  void apfFEA::Adapt()
  {
    std::cerr << "APF Adaption loop needs to be added!\n";
    exit(1);
  }
}  // namespace amsi
