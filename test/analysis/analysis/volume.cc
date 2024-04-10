#include "amsiAnalysis.h"
#include "apfMeasure.h"
#include <cstdlib>
int main(int argc, char * argv[])
{
  int result = 0;
  //std::string lic = std::getenv("SIM_LICENSE_FILE");
  //amsi::useSimmetrix(lic);
  //amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  //pGModel mdl = GM_load(argv[1],0,NULL);
  //pParMesh msh = PM_load(argv[2],mdl,NULL);
  //int tg = atoi(argv[3]);
  //apf::Mesh * apf_msh = apf::createMesh(msh);
  //apf::Field * u =  apf::createLagrangeField(apf_msh,"u",apf::VECTOR,1);
  //apf::zeroField(u);
  //pGEntity rgn = (pGRegion)GM_entityByTag(mdl,3,tg);
  //double vol_surf = 0.0;
  //vol_surf = amsi::measureDisplacedModelEntity_greens(reinterpret_cast<apf::ModelEntity*>(rgn),u);
  //double vol_rgn = amsi::measureDisplacedModelEntity(reinterpret_cast<apf::ModelEntity*>(rgn),u);
  //int rnk = -1;
  //MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
  //if(rnk == 0)
  //  std::cout << vol_surf << " " << vol_rgn << std::endl;
  //amsi::freeAnalysis();
  return result;
}
