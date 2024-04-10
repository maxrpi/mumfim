#ifndef AMSI_APF_MEASURE_H_
#define AMSI_APF_MEASURE_H_
#include <apf.h>
#include <apfGeometry.h>
#include <apfMesh.h>
namespace amsi
{
  double distance(const apf::Vector3 & pt_a, const apf::Vector3 & pt_b);
  double triangleArea(const apf::Vector3 & pt_a, const apf::Vector3 & pt_b, const apf::Vector3 & pt_c);
  void faceNormal(const apf::Vector3 & pt_a, const apf::Vector3 & pt_b, const apf::Vector3 & pt_c, apf::Vector3 & n);
  int side(apf::ModelEntity * mdl_ent, apf::Mesh * msh, apf::MeshEntity * fc);
  double measureModelEntity(apf::ModelEntity * mdl_ent, apf::Mesh * msh);
  double measureDisplacedModelEntity(apf::ModelEntity * mdl_ent, apf::Field * u);
  /**
   * @brief Measure the mesh entity w.r.t the entity's dimension, length, area, or volume
   * @param ent The MeshEntity to measure
   * @param u The displacement Field to apply to ent
   * @note Only implemented for regions which gives volume at the moment
   */
  double measureDisplacedMeshEntity(apf::MeshEntity * ent, apf::Field * u);
  double measureDisplacedModelEntity_greens(apf::ModelEntity * mdl_ent, apf::Field * u);
  template <typename I>
    double measureDisplacedModelEntities(I b, I e, apf::Field * u);
}
#include "apfMeasure_impl.h"
#endif
