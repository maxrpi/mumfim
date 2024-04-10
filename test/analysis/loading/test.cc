#include "amsiAnalysis.h"
#include <cassert>
#include <cstdlib>
int main(int argc, char * argv[])
{
  assert(argc == 3);
  int result = 0;
  amsi::initAnalysis(argc,argv, MPI_COMM_WORLD);
  amsi::freeAnalysis();
  return result;
}
