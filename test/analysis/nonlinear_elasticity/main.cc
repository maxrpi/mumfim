#include "amsiAnalysis.h"
#include "NonLinElasticity.h"
#include "Solvers.h"
#include <mpi.h>
#include <cassert>
#include <iostream>
int main (int argc, char ** argv)
{
  int result = 0;
  amsi::useSimmetrix("/net/common/meshSim/license/license.txt");
  amsi::usePetsc("petsc_options");
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  Sim_logOn("sim.log");
  pGModel mdl = GM_load(argv[1],0,NULL);
  pParMesh msh = PM_load(argv[2],mdl,NULL);
  std::vector<pACase> css;
  amsi::getTypeCases(SModel_attManager(mdl),"analysis",std::back_inserter(css));
  amsi::initCase(mdl,css[0]);
  pACase pd = (pACase)AttNode_childByType(
      (pANode)css[0], amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
  pACase ss = (pACase)AttNode_childByType(
      (pANode)css[0], amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
  amsi::LAS * las = new amsi::PetscLAS(0,0);
  apf::Mesh * apf_msh = apf::createMesh(msh);
  std::vector<apf::Field*> flds;
  amsi::buildFieldsFromSim(pd,apf_msh,std::back_inserter(flds));
  std::vector<apf::Numbering*> nms;
  amsi::buildNumberingsFromSim(pd,flds.begin(),flds.end(),std::back_inserter(nms));;
  std::vector<amsi::SimBCQuery*> bcs[amsi::BCType::num_bc_types]; // can't init on the BGQ
  amsi::buildBCQueriesFromSim(pd,amsi::BCType::dirichlet,std::back_inserter(bcs[amsi::BCType::dirichlet]));
  std::vector<pANode> bc_nds;
  amsi::getTypeNodes((pANode)pd,amsi::getBCTypeString(amsi::BCType::dirichlet),std::back_inserter(bc_nds));
  for(auto bc_nd = bc_nds.begin(); bc_nd != bc_nds.end(); ++bc_nd)
  {
    amsi::describeNode(*bc_nd);
    std::vector<pANode> fld_nm;
    amsi::getTypeNodes(*bc_nd,"field name",std::back_inserter(fld_nm));
    amsi::describeNode(fld_nm[0]);
  }
  //amsi::buildSimBCQueries(pd,amsi::DIRICHLET,/*amsi::UNITLESS,amsi::NUM_FIELD_UNITS*/,std::back_inserter(bcs[amsi::DIRICHLET]));
  //amsi::buildSimBCQueries(pd,amsi::NEUMANN,/*amsi::CUSTOM,amsi::NUM_NEUMANN_TYPES*/,std::back_inserter(bcs[amsi::NEUMANN]));
  //after building queries, associate them with fields, since we need the geometric entities they lie on anyway, but we also have to relate the SimBCQueries BACK to the attributes taged on the model...
  amsi::NonLinElasticity iso_non(mdl,msh,pd,ss);
  double nrm = 0.0;
  amsi::NewtonSolver(&iso_non,las,30,1e-4,1.0,nrm);
  apf::writeVtkFiles("isotropic_nonlinear_elastic_result",iso_non.getMesh());
  delete las;
  amsi::freeCase(css[0]);
  Sim_logOff();
  amsi::freeAnalysis();
  return result;
}
