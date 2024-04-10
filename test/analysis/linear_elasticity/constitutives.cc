#include <amsiLinearElasticConstitutive.h>
#include <mpi.h>
int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  apf::DynamicMatrix C = amsi::isotropicLinearElasticityTensor(10000,0.3);
  //apf::Field * fld = apf::createUserField(NULL,"test_field",apf::SCALAR,
  //amsi::LinearElasticIntegrator * itgr = new amsi::LinearElasticIntegrator(fld,1);
  //delete itgr;
  MPI_Finalize();
  return 0;
}
