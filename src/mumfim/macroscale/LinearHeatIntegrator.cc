#include "LinearHeatIntegrator.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <cstring> //memset
#include <cmath>
namespace mumfim
{

  apf::Matrix3x3 * LinearHeatIntegrator::getIsotropicKappa(double k)
  {
    apf::Matrix3x3 * K = new apf::Matrix3x3(
      k   ,  0.0 ,  0.0 ,
      0.0 ,    k ,  0.0 ,
      0.0 ,  0.0 ,    k
    );
    return K;
  }

  apf::Matrix3x3 * LinearHeatIntegrator::getFullKappa(double xx, double yy, double zz,
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

  LinearHeatIntegrator::LinearHeatIntegrator(apf::Field * temperature, apf::Matrix3x3 * kappa)
    : amsi::ElementalSystem(temperature,1)
    , T_(temperature)
    , D(apf::fromMatrix(*kappa))
  {
    ;
  }

  LinearHeatIntegrator::LinearHeatIntegrator(apf::Field * temperature, apf::Field * kappa)
    : amsi::ElementalSystem(temperature,1)
    , T_(temperature)
    , K_(kappa)
  {
    ;
  }

  void LinearHeatIntegrator::inElement(apf::MeshElement *me)
  {
    // Populate kappa_e (elemental kappa) with IP value
    // from K_ field for this integrator interation
    //apf::Matrix3x3 K[1];
    //K_->getNodeValue(apf::getMeshEntity(me), 0, K);
    //D = apf::fromMatrix<3,3>(K[0]);
    ;
  }

  void LinearHeatIntegrator::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    // Using TJHughes, pp 68-70, add nonhomogenous part later
    //apf::NewArray<double> Na; apf::getShapeValues(e, p, Na);
    apf::NewArray<apf::Vector3> Ba; apf::getShapeGrads(e, p, Ba);

    apf::DynamicMatrix B(3, nedofs);
    apf::DynamicMatrix BT(nedofs, 3);
    for(int i = 0; i < nedofs; i++)
    {
      B(0, i) = BT(i, 0) = Ba[i][0];
      B(1, i) = BT(i, 1) = Ba[i][1];
      B(2, i) = BT(i, 2) = Ba[i][2];
    }

    apf::DynamicMatrix BT_D_B(nedofs,nedofs);
    apf::DynamicMatrix D_B(3,nedofs); // intermediate product (the flux!)
    apf::multiply(D, B, D_B);
    apf::multiply(BT, D_B, BT_D_B);
    // numerical integration
    BT_D_B *= w * dV;  // BT_D_B now elemental integral contribution for p (ke_p)
    Ke += BT_D_B;  // Accumulate over integration points
  }
}