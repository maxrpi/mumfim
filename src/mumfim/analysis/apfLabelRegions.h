#ifndef APF_SIM_WRAPPER_H_
#define APF_SIM_WRAPPER_H_
#include <apf.h>
#include <apfMesh.h>
namespace amsi
{
  /*
   * \brief applies a integer tag to each region
   */
  void applyUniqueRegionTags(apf::Mesh* mesh);
}
#endif
