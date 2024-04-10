#include "test.h"
#include <apfMeshGen.h>
int main(int argc, char * argv[])
{
  int fld = 0;
  apf::Mesh2 * cbe = amsi::fiveTetCube();
  apf::writeOneVtkFile("cube",cbe);
  return fld;
}
