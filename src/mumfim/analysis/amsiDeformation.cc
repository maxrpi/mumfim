#include "amsiDeformation.h"
namespace amsi {
  const apf::Matrix3x3 eye3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  /*
   * \param[in] e - mesh element with the displacements as the underlying field
   * \param[in] p - local coordinate where the deformation gradient is
   * \param[out] F - computed deformation gradient
   */
  void deformationGradient(apf::Element* e, const apf::Vector3& p,
                           apf::Matrix3x3& F)
  {
    apf::NewArray<apf::Vector3> u;
    apf::getVectorNodes(e, u);
    apf::NewArray<apf::Vector3> Ni_j;
    apf::getShapeGrads(e, p, Ni_j);
    F = eye3x3;
    int nen = apf::countNodes(e);
    // get vector nodes always has 3 components
    int cmp = 3;
    for (int ii = 0; ii < cmp; ++ii)
      for (int jj = 0; jj < cmp; ++jj)
        for (int kk = 0; kk < nen; ++kk) F[ii][jj] += u[kk][ii] * Ni_j[kk][jj];
  }
  /*
  void deformationGradientDerivative(apf::Element* e, const apf::Vector3& p,
                           apf::NewArray<apf::Matrix3x3>& dFdxi)
  {
    //number of nodes can be found by shape->coundNodes;
    //number of components can be found by field->coundComponents;
    apf::NewArray<apf::Vector3> Ni_j;
    apf::getShapeGrads(e, p, Ni_j);
    // for each node: dFdxi = Fi*Ni_j
  }
  */
}  // namespace amsi
