#include "apfMatrixUtil.h"
#include <cassert>
namespace amsi
{
  apf::Matrix3x3 symmetricPart(apf::Matrix3x3 & m)
  {
    apf::Matrix3x3 rslt;
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        rslt[ii][jj] = 0.5 * (m[ii][jj] + m[jj][ii]);
    return rslt;
  }
  apf::Matrix3x3 antisymmetricPart(apf::Matrix3x3 & m)
  {
    apf::Matrix3x3 rslt;
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        rslt[ii][jj] = 0.5 * (m[ii][jj] - m[jj][ii]);
    return rslt;
  }
  void voigtVec2Mat(int dim, const double * vec, apf::Matrix3x3 & mat)
  {
    if(dim == 2)
    {
      mat[0][0] = vec[0];
      mat[1][1] = vec[1];
      mat[0][1] = vec[2];
    }
    else if(dim == 3)
    {
      mat[0][0] = vec[0];
      mat[1][1] = vec[1];
      mat[2][2] = vec[2];
      mat[1][2] = vec[3];
      mat[0][2] = vec[4];
      mat[0][1] = vec[5];
    }
  }
  void mat2VoigtVec(int dim, const apf::Matrix3x3 & mat, double * vec)
  {
    if(dim == 2)
    {
      vec[0] = mat[0][0];
      vec[1] = mat[1][1];
      vec[2] = mat[0][1];
    }
    else if(dim == 3)
    {
      vec[0] = mat[0][0];
      vec[1] = mat[1][1];
      vec[2] = mat[2][2];
      vec[3] = mat[1][2];
      vec[4] = mat[0][2];
      vec[5] = mat[0][1];
    }
  }
  int mat2VoigtIdx(int dim, int rr, int cc)
  {
    assert(dim == 2 || dim == 3);
    static const int mp2d[2][2] = {{0,2},{2,1}};
    static const int mp3d[3][3] = {{0,5,4},{5,1,3},{4,3,2}};
    if(dim == 2)
      return mp2d[rr][cc];
    else if (dim == 3)
      return mp3d[rr][cc];
    return -1;
  }
  int voigt2MatRow(int dim, int vgt)
  {
    assert(dim == 2 || dim == 3);
    static const int mp2d[] = {0,1,0};
    static const int mp3d[] = {0,1,2,1,0,0};
    if(dim == 2)
      return mp2d[vgt];
    else if(dim == 3)
      return mp3d[vgt];
    return -1;
  }
  int voigt2MatCol(int dim, int vgt)
  {
    assert(dim == 2 || dim == 3);
    static const int mp2d[] = {0,1,0};
    static const int mp3d[] = {0,1,2,2,2,1};
    if(dim == 2)
      return mp2d[vgt];
    else if (dim == 3)
      return mp3d[vgt];
    return -1;
  }
  void mat2Array(apf::DynamicMatrix & mat, double * ptr)
  {
    int rcnt = mat.getRows();
    int ccnt = mat.getColumns();
    for(int rr = 0; rr < rcnt; ++rr)
      for(int cc = 0; cc < ccnt; ++cc)
        ptr[rr*ccnt+cc] = mat(rr,cc);
  }
}
