#include "LinearHeatIntegrator.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <cstring> //memset
#include <cmath>
#include <cassert>
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

  LinearHeatIntegrator::LinearHeatIntegrator(apf::Field * temperature, apf::Numbering* numbering, apf::Matrix3x3 * kappa)
    : amsi::ElementalSystem(temperature,numbering, 1)
    , T_(temperature)
    , D(apf::fromMatrix(*kappa))
  {
    ;
  }

  LinearHeatIntegrator::LinearHeatIntegrator(apf::Field * temperature,apf::Numbering* numbering, apf::Field * kappa)
    : amsi::ElementalSystem(temperature,numbering, 1)
    , T_(temperature)
    , K_(kappa)
  {
    ;
  }

  /*
  void LinearHeatIntegrator::inElement(apf::MeshElement *me)
  {
    e = apf::createElement(T_, me);
    nenodes = apf::countNodes(e);
    nedofs = nenodes;  // Assumes scalar values at nodes
    Ke.setSize(nenodes, nenodes);
    Ke.zero();
    fe.setSize(nenodes);
    fe.zero();
  }
  */

  void LinearHeatIntegrator::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    // Using TJHughes, pp 68-70, add nonhomogenous part later
    // fe holds the residual vector, because of the incremental formulation
    //apf::NewArray<double> Na; apf::getShapeValues(e, p, Na);
    apf::NewArray<apf::Vector3> Ba; apf::getShapeGrads(e, p, Ba);

    apf::DynamicMatrix B(3, nenodes);
    apf::DynamicMatrix BT(nenodes, 3);
    apf::DynamicVector Theta(nenodes);
    assert(dV>0);
    for(int i = 0; i < nenodes; i++)
    {
      B(0, i) = BT(i, 0) = Ba[i][0];
      B(1, i) = BT(i, 1) = Ba[i][1];
      B(2, i) = BT(i, 2) = Ba[i][2];
      Theta(i) = field_values_[i];
    }

    apf::DynamicMatrix BT_D_B(nenodes,nenodes);
    apf::DynamicMatrix D_B(3,nenodes); // intermediate product (the flux!)
    apf::multiply(D, B, D_B);
    apf::multiply(BT, D_B, BT_D_B);
    // numerical integration
    BT_D_B *= w * dV;  // BT_D_B now elemental integral contribution for p (ke_p)
    Ke += BT_D_B;  // Accumulate over integration points
    // fe holds the residual vector, because of the incremental formulation
    // In the case of heat sources the residual will be Ke*theta - rhs
    apf::DynamicVector fe_ip(nedofs);
    apf::multiply(Ke, Theta, fe_ip); 
    fe -=  fe_ip;
  }
}