#include "amsiAnalysis.h"
#include "simAttributes.h"
#include "amsiElasticityFEA.h"
#include "Solvers.h"
#include <mpi.h>
#include <cassert>
#include <iostream>
int main (int argc, char ** argv)
{
  assert(argc == 3);
  int result = 0;
  std::string lic = std::getenv("SIM_LICENSE_FILE");
  amsi::useSimmetrix(std::string(lic));
  amsi::usePetsc(std::string("petsc_options"));
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  {
    Sim_logOn("sim.log");
    pGModel mdl = GM_load(argv[1],0,NULL);
    pParMesh msh = PM_load(argv[2],mdl,NULL);
    std::vector<pACase> css;
    amsi::getTypeCases(SModel_attManager(mdl),"analysis",std::back_inserter(css));
    amsi::initCase(mdl,css[0]);
    pACase pd = (pACase)AttNode_childByType(
        (pANode)css[0],
        amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
    pACase ss = (pACase)AttNode_childByType(
        (pANode)css[0],
        amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
    amsi::PetscLAS las(0,0);
    amsi::ElasticityFEA iso_lin(mdl,msh,pd,ss,AMSI_COMM_SCALE);
    amsi::LinearSolver(&iso_lin,&las);
    apf::writeVtkFiles("isotropic_linear_elastic_result",
                       iso_lin.getMesh());
    amsi::freeCase(css[0]);
    Sim_logOff();
  } // destroys all stack-allocated objects before deinitialization of libraries
  amsi::freeAnalysis();
  return result;
}
