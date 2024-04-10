#ifndef AMSI_DEFORMATION_H_
#define AMSI_DEFORMATION_H_
#include <apf.h>
namespace amsi
{
  void deformationGradient(apf::Element * e, const apf::Vector3 & p, apf::Matrix3x3 & F);
  // derivative of deformation gradient at point w.r.t. each node
  void deformationGradientDerivative(apf::Element * e, const apf::Vector3 & p,
                                     apf::NewArray<apf::Matrix3x3> & dFdxi);
}
#endif
