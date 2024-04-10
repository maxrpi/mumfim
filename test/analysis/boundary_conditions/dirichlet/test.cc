#include "test.h"
#include "amsiAnalysis.h"
#include "simAnalysis.h"
#include "simAttributes.h"
#include "simBoundaryConditions.h"
#include "apf.h"
#include "apfSIM.h"
#include "apfBoundaryConditions.h"
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
  pGModel mdl = GM_load(argv[1],0,NULL);
  pParMesh sm_msh = PM_load(argv[2],mdl,NULL);
  pMesh prt = PM_mesh(sm_msh,0);
  // get all analysis cases
  std::vector<pACase> css;
  amsi::getTypeCases(SModel_attManager(mdl),"analysis",std::back_inserter(css));
  // only run the first attribute case
  amsi::initCase(mdl,css[0]);
  pACase pd = (pACase)AttNode_childByType((pANode)css[0],"problem definition");
  int dsp = amsi::FieldUnit::displacement;
  apf::Mesh * msh =  apf::createMesh(sm_msh);
  apf::Field * u = apf::createLagrangeField(msh,"displacement",apf::VECTOR,1);
  apf::Numbering * nm = apf::createNumbering(u);
  // wrap below here in function ?
  std::vector<amsi::SimBCQuery*> dir_qrys;
  amsi::buildSimBCQueries(pd,amsi::BCType::dirichlet,&dsp,(&dsp)+1,std::back_inserter(dir_qrys));
  int fxd = amsi::applySimDirichletBCs(nm,prt,dir_qrys.begin(),dir_qrys.end(),0.0);
  failed += test_neq("Fixed dofs",0,fxd);
  amsi::freeCase(css[0]);
  amsi::freeAnalysis();
  return failed;
}
