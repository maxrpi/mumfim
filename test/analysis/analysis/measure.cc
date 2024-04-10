#include <apfMeasure.h>
#include <apfMeshUtil.h>
#include <amsiAnalysis.h>
#include <cmath>
inline bool close(double a, double b, double eps = 1e-8)
{
  return abs(a) - abs(b) < eps;
}
int main(int argc, char * argv[])
{
  int result = 0;
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  /*
  apf::Vector3 pt(0.0,0.0,0.0);
  apf::Mesh2 * pt_msh = amsi::makeNullMdlEmptyMesh();
  pt_msh->createVertex(NULL,pt,apf::Vector3(0.0,0.0,0.0));
  pt_msh->acceptChanges();
  apf::Field * pt_msh_u = apf::createLagrangeField(pt_msh,"u",apf::VECTOR,1);
  apf::zeroField(pt_msh_u);
  apf::MeshEntity * pt_ent = amsi::getFirstMeshEntity(pt_msh,0);
  double pt_msr = amsi::measureDisplacedMeshEntity(pt_ent,pt_msh_u);
  result += pt_msr == 0.0 ? 0 : 1;
  */
  apf::Vector3 ln[] = {apf::Vector3(0.0,0.0,0.0), apf::Vector3(1.0,0.0,0.0)};
  apf::Mesh * ln_msh = amsi::makeSingleEntityMesh(apf::Mesh::EDGE,&ln[0]);
  apf::Field * ln_msh_u = apf::createLagrangeField(ln_msh,"u",apf::VECTOR,1);
  apf::zeroField(ln_msh_u);
  apf::MeshEntity * ln_ent = amsi::getFirstMeshEntity(ln_msh,1);
  double ln_msr = amsi::measureDisplacedMeshEntity(ln_ent,ln_msh_u);
  result += close(ln_msr,1.0) ? 0 : 1;
  apf::Vector3 tri[] = {apf::Vector3(0.0,0.0,0.0), apf::Vector3(1.0,0.0,0.0), apf::Vector3(0.0,1.0,0.0)};
  apf::Mesh * tri_msh = amsi::makeSingleEntityMesh(apf::Mesh::TRIANGLE,&tri[0]);
  apf::Field * tri_msh_u = apf::createLagrangeField(tri_msh,"u",apf::VECTOR,1);
  apf::zeroField(tri_msh_u);
  apf::MeshEntity * tri_ent = amsi::getFirstMeshEntity(tri_msh,2);
  double tri_msr = amsi::measureDisplacedMeshEntity(tri_ent,tri_msh_u);
  result += close(tri_msr,0.5) ? 0 : 1;
  apf::Vector3 tet[] = {apf::Vector3(0.0,0.0,0.0), apf::Vector3(1.0,0.0,0.0), apf::Vector3(0.0,0.0,1.0), apf::Vector3(0.0,1.0,0.0)};
  apf::Mesh * tet_msh = amsi::makeSingleEntityMesh(apf::Mesh::TET,&tet[0]);
  apf::Field * tet_msh_u = apf::createLagrangeField(tet_msh,"u",apf::VECTOR,1);
  apf::zeroField(tet_msh_u);
  apf::MeshEntity * tet_ent = amsi::getFirstMeshEntity(tet_msh,3);
  double tet_msr = amsi::measureDisplacedMeshEntity(tet_ent,tet_msh_u);
  result += close(tet_msr,0.25) ? 0 : 1;
  amsi::freeAnalysis();
  return result;
}
