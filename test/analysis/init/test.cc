#include "test.h"
#include "amsiAnalysis.h"
int main(int argc, char ** argv)
{
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  amsi::freeAnalysis();
  return 0;
}
