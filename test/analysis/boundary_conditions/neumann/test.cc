#include "test.h"
#include "amsiAnalysis.h"
#include "amsiDummyLAS.h"
#include "simAnalysis.h"
#include "simAttributes.h"
#include "simBoundaryConditions.h"
#include "apfBoundaryConditions.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfSIM.h>
#include <mpi.h>
#include <cassert>
#include <fstream>
int main(int argc, char ** argv)
{
  int failed = 0;
  assert(argc == 3);
  std::string lic = std::getenv("SIM_LICENSE_FILE");
  amsi::useSimmetrix(lic);
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  Sim_logOn("simlog");
  pGModel mdl = GM_load(argv[1],0,NULL);
  pParMesh sm_msh = PM_load(argv[2],mdl,NULL);
  pMesh prt = PM_mesh(sm_msh,0);
  (void)prt;
  // get all analysis cases
  std::vector<pACase> css;
  amsi::getTypeCases(SModel_attManager(mdl),"analysis",std::back_inserter(css));
  // only run the first attribute case
  amsi::initCase(mdl,css[0]);
  pACase pd = (pACase)AttNode_childByType((pANode)css[0],"problem definition");
  (void)pd;
  apf::Mesh * msh =  apf::createMesh(sm_msh);
  apf::Field * u = apf::createLagrangeField(msh,"displacement",apf::VECTOR,1);
  apf::Numbering * nm = apf::createNumbering(u);
  apf::NaiveOrder(nm);
  /*
  amsi::DummyLAS las(dofs);
  int tps[] = {amsi::NeuBCType::traction,amsi::NeuBCType::pressure};
  std::vector<amsi::SimBCQuery*> neu_qrys;
  amsi::buildSimBCQueries(pd,amsi::BCType::neumann,&tps[0],(&tps[0])+1,std::back_inserter(neu_qrys));
  amsi::applySimNeumannBCs(&las,nm,prt,neu_qrys.begin(),neu_qrys.end(),1.0);
  double nrm = 0.0;
  las.GetVectorNorm(nrm);
  failed += test_neq("Force vector norm",0.0,nrm);
  */
  amsi::freeCase(css[0]);
  Sim_logOff();
  amsi::freeAnalysis();
  return failed;
}
