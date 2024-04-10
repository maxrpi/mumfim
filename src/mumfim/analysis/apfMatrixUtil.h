#ifndef APF_MATRIX_UTIL_H_
#define APF_MATRIX_UTIL_H_
#include <apfMatrix.h>
#include <apfDynamicMatrix.h>
namespace amsi
{
  apf::Matrix3x3 symmetricPart(apf::Matrix3x3 & m);
  apf::Matrix3x3 antisymmetricPart(apf::Matrix3x3 & m);
  void voigtVec2Mat(int dim, const double * vec, apf::Matrix3x3 & mat);
  void mat2VoigtVec(int dim, const apf::Matrix3x3 & mat, double * vec);
  int mat2VoigtIdx(int dim, int rr, int cc);
  int voigt2MatRow(int dim, int vgt);
  int voigt2MatCol(int dim, int vgt);
  void mat2Array(apf::DynamicMatrix & mat, double * ptr);
}
#endif
