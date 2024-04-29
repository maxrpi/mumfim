#include "amsiFEA.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <vector>
namespace amsi {
  FEAStep::FEAStep(const ModelDefinition& pd,
           const ModelDefinition& ss,
           const ModelDefinition& out, std::vector<DirichletBCEntry> dbc,
           std::vector<NeumannBCEntry> nbc, const std::string& analysis_name,
           MPI_Comm cm)
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
  {
  }

  FEAStep::FEAStep(const mt::CategoryNode& analysis_case,
           std::vector<DirichletBCEntry> dbc, std::vector<NeumannBCEntry> nbc,
           const std::string& analysis_name, MPI_Comm cm)
      : constraint_dofs(0)
      , local_dof_count(0)
      , first_local_dof(0)
      , global_dof_count(0)
      , fixed_dofs(0)
      , numbered(false)
      , T(0.0)
      , analysis_comm(cm)
      , dirichlet_bcs(std::move(dbc))
      , neumann_bcs(std::move(nbc))
  {
    const mt::CategoryNode* pd;
    const mt::CategoryNode* ss;
    const mt::CategoryNode* out;
    const mt::CategoryNode* tmp;
    pd = analysis_case.FindCategoryNodeByType("problem definition");
    out = analysis_case.FindCategoryNodeByType("output");
    ss = analysis_case.FindCategoryNodeByType("solution strategy");
    if (pd == nullptr || ss == nullptr || out == nullptr) {
      std::cerr << "analysis case should have \"problem "
                   "definition\",\"output\", and \"solution strategy\".\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    tmp = pd->FindCategoryNodeByType(analysis_name);
    pd = (tmp == nullptr) ? pd : tmp;
    tmp = ss->FindCategoryNodeByType(analysis_name);
    ss = (tmp == nullptr) ? ss : tmp;
    tmp = out->FindCategoryNodeByType(analysis_name);
    out = (tmp == nullptr) ? out : tmp;
    problem_definition = ModelDefinition{
        .associated = mt::AssociateModel(pd),
        .unassociated = std::make_shared<mt::CategoryNode>(*pd)};
    solution_strategy = ModelDefinition{
        .associated = mt::AssociateModel(ss),
        .unassociated = std::make_shared<mt::CategoryNode>(*ss)};
    output = ModelDefinition{
        .associated = mt::AssociateModel(out),
        .unassociated = std::make_shared<mt::CategoryNode>(*out)};
    assert(solution_strategy.associated->GetNullGeometry() != nullptr);
  }
  void FEAStep::setSimulationTime(double t)
  {
    T = t;
    if (!PCU_Comm_Self())
      std::cout << "Simulation time updated: " << T << std::endl;
  }
  void FEAStep::GetDOFInfo(int& global, int& local, int& offset)
  {
    global = global_dof_count;
    local = local_dof_count;
    offset = first_local_dof;
  }
  void assembleMatrix(LAS* las, int rw_cnt, int* rw_nms, int cl_cnt,
                      int* cl_nms, double* Ke)
  {
    las->AddToMatrix(rw_cnt, rw_nms, cl_cnt, cl_nms, &Ke[0]);
  }
  void assembleVector(LAS* las, int rw_cnt, int* rw_nms, double* fe)
  {
    las->AddToVector(rw_cnt, rw_nms, &fe[0]);
  }
}  // namespace amsi
