#include "apfWrapper.h"
namespace amsi
{
  void getScalarDofValues(apf::Element * elmt, double * vls)
  {
    int nenodes = apf::countNodes(elmt);
    apf::NewArray<double> vs;
    apf::getScalarNodes(elmt,vs);
    for(int ii = 0; ii < nenodes; ii++)
      vls[ii] = vs[ii];
  }
  void getVectorDofValues(apf::Element * elmt, double * vls)
  {
    int nenodes = apf::countNodes(elmt);
    apf::NewArray<apf::Vector3> vs;
    apf::getVectorNodes(elmt,vs);
    for(int ii = 0; ii < nenodes; ii++)
      for(int jj = 0; jj < 3; jj++)
        vls[ii*3 + jj] = vs[ii][jj];
  }
  void getMatrixDofValues(apf::Element * elmt, double * vls)
  {
    int nenodes = apf::countNodes(elmt);
    apf::NewArray<apf::Matrix3x3> vs;
    apf::getMatrixNodes(elmt,vs);
    for(int ii = 0; ii < nenodes; ii++)
      for(int jj = 0; jj < 3; jj++)
        for(int kk = 0; kk < 3; kk++)
          vls[ii*9 + ii*3 + kk] = vs[ii][jj][kk];
  }
  void getDofValues(apf::Field * fld, apf::Element * elmt, double * vls)
  {
    int tp = apf::getValueType(fld);
    switch(tp)
    {
    case apf::SCALAR:
      getScalarDofValues(elmt,vls);
      break;
    case apf::VECTOR:
      getVectorDofValues(elmt,vls);
      break;
    case apf::MATRIX:
      getMatrixDofValues(elmt,vls);
      break;
    }
  }
}
