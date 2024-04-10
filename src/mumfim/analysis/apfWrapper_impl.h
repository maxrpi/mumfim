#include <amsiMPI.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <iostream>
namespace amsi
{
  template <class I>
    void getAvgScalarFieldValue(apf::Field * fld, I bgn, I nd, double * avg)
  {
    *avg = 0.0;
    long cnt = 0;
    apf::FieldShape * fs = apf::getShape(fld);
    apf::Mesh * msh = apf::getMesh(fld);
    for(I it = bgn; it != nd; ++it)
    {
      apf::MeshEntity * ent = reinterpret_cast<apf::MeshEntity*>(*it);
      int nd_cnt = fs->countNodesOn(msh->getType(ent));
      cnt += nd_cnt;
      for(int ii = 0; ii < nd_cnt; ii++)
        *avg += apf::getScalar(fld,ent,ii);
    }
    *avg = comm_avg(*avg,cnt);
  }
  template <class I>
    void getAvgVectorFieldValue(apf::Field * fld, I bgn, I nd, double * avg)
  {
    for(int ii = 0; ii < 3; ii++)
      avg[ii] = 0.0;
    long cnt = 0;
    apf::FieldShape * fs = apf::getShape(fld);
    apf::Mesh * msh = apf::getMesh(fld);
    for(I it = bgn; it != nd; ++it)
    {
      apf::MeshEntity * ent = reinterpret_cast<apf::MeshEntity*>(*it);
      int nd_cnt = fs->countNodesOn(msh->getType(ent));
      cnt += nd_cnt;
      apf::Vector3 val;
      for(int ii = 0; ii < nd_cnt; ii++)
      {
        apf::getVector(fld,ent,ii,val);
        for(int jj = 0; jj < 3; jj++)
          avg[jj] += val[jj];
      }
    }
    for(int ii = 0; ii < 3; ii++)
      avg[ii] = comm_avg(avg[ii],cnt);
  }
  template <class I>
    void getAvgMatrixFieldValue(apf::Field * fld, I bgn, I nd, double * avg)
  {
    for(int ii = 0; ii < 9; ii++)
      avg[ii] = 0.0;
    long cnt = 0;
    apf::FieldShape * fs = apf::getShape(fld);
    apf::Mesh * msh = apf::getMesh(fld);
    for(I it = bgn; it != nd; ++it)
    {
      apf::MeshEntity * ent = reinterpret_cast<apf::MeshEntity*>(*it);
      int nd_cnt = fs->countNodesOn(msh->getType(ent));
      cnt += nd_cnt;
      apf::Matrix3x3 val;
      for(int ii = 0; ii < nd_cnt; ii++)
      {
        apf::getMatrix(fld,ent,ii,val);
        for(int jj = 0; jj < 9; jj++)
          avg[jj] += val[jj / 3][jj % 3];
      }
    }
    for(int ii = 0 ; ii < 9; ii++)
      avg[ii] = comm_avg(avg[ii],cnt);
  }
  template <class I>
    void getAvgFieldValue(apf::Field * fld, I bgn, I nd, double * avg)
  {
    int tp = apf::getValueType(fld);
    switch(tp)
    {
    case apf::SCALAR:
      getAvgScalarFieldValue(fld,bgn,nd,avg);
      break;
    case apf::VECTOR:
      getAvgVectorFieldValue(fld,bgn,nd,avg);
      break;
    case apf::MATRIX:
      getAvgMatrixFieldValue(fld,bgn,nd,avg);
      break;
    default:
      std::cerr << "ERROR: unrecognized field type " << tp << std::endl;
    }
  }
}
