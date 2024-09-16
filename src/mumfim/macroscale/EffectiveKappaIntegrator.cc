#include "EffectiveKappaIntegrator.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <cstring> //memset
#include <cmath>
#include <cassert>
namespace mumfim
{

  apf::Matrix3x3 * EffectiveKappaIntegrator::getIsotropicKappa(double k)
  {
    apf::Matrix3x3 * K = new apf::Matrix3x3(
      k   ,  0.0 ,  0.0 ,
      0.0 ,    k ,  0.0 ,
      0.0 ,  0.0 ,    k
    );
    return K;
  }

  apf::Matrix3x3 * EffectiveKappaIntegrator::getFullKappa(double xx, double yy, double zz,
    double xy, double yz, double zx)
  {
    //Onsager reciprocity -> symmetrical
    apf::Matrix3x3 * K = new apf::Matrix3x3(
      xx , xy , zx ,
      xy , yy , yz ,
      zx , yz , zz
    );
    return K;
  }

  EffectiveKappaIntegrator::EffectiveKappaIntegrator(apf::Field *dummy, apf::Numbering* numbering, apf::Matrix3x3 * kappa,
    std::map<std::pair<int, int>, double *> &k_map_)
    : amsi::ElementalSystem(dummy,numbering, 1)
    , D(apf::fromMatrix(*kappa))
    , k_map(k_map_)
  {
    ;
  }


  void EffectiveKappaIntegrator::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    // Using TJHughes, pp 68-70, add nonhomogenous part later
    // fe holds the residual vector, because of the incremental formulation
    //apf::NewArray<double> Na; apf::getShapeValues(e, p, Na);
    apf::NewArray<apf::Vector3> Ba; apf::getShapeGrads(e, p, Ba);

    apf::DynamicMatrix B(3, nenodes);
    apf::DynamicMatrix BT(nenodes, 3);
    assert(dV>0);
    for(int i = 0; i < nenodes; i++)
    {
      B(0, i) = BT(i, 0) = Ba[i][0];
      B(1, i) = BT(i, 1) = Ba[i][1];
      B(2, i) = BT(i, 2) = Ba[i][2];
    }

    apf::DynamicMatrix BT_D_B(nenodes,nenodes);
    apf::DynamicMatrix D_B(3,nenodes);
    apf::multiply(D, B, D_B);
    apf::multiply(BT, D_B, BT_D_B);

    BT_D_B *= w * dV;  

    for (int a=0; a < nenodes; a++){
      int global_na = field_numbers_[a];
      for (int b=0; b < nenodes; b++){
        int global_nb = field_numbers_[b];
        std::pair<int, int> key(global_na, global_nb);
        *k_map[key] += BT_D_B(a,b);
      }
    }
  }
}