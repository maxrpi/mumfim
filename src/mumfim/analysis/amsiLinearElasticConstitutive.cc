#include "amsiLinearElasticConstitutive.h"
#include <apfMesh.h>
#include <apfShape.h>
#include <cstring> //memset
#include <cmath>
namespace amsi
{
  apf::DynamicMatrix isotropicLinearElasticityTensor(double E, double v)
  {
    apf::DynamicMatrix result(6,6);
    double lambda = ( v * E ) / ( ( 1 + v ) * ( 1 - 2 * v ) );
    double mu = E / ( 2 * ( 1 + v ) );
    result(0,0) = lambda + (2 * mu);
    result(0,1) = lambda;
    result(0,2) = lambda;
    result(0,3) = result(0,4) = result(0,5) = 0.0;
    result(1,0) = lambda;
    result(1,1) = lambda + (2 * mu);
    result(1,2) = lambda;
    result(1,3) = result(1,4) = result(1,5) = 0.0;
    result(2,0) = lambda;
    result(2,1) = lambda;
    result(2,2) = lambda + (2 * mu);
    result(2,3) = result(2,4) = result(2,5) = 0.0;
    result(3,0) = result(3,1) = result(3,2) = result(3,4) = result(3,5) = 0.0;
    result(3,3) = mu;
    result(4,0) = result(4,1) = result(4,2) = result(4,3) = result(4,5) = 0.0;
    result(4,4) = mu;
    result(5,0) = result(5,1) = result(5,2) = result(5,3) = result(5,4) = 0.0;
    result(5,5) = mu;
    return result;
  }
  LinearElasticIntegrator::LinearElasticIntegrator(apf::Field * field, apf::Numbering* numbering, int o, double E, double v)
    : ElementalSystem(field, numbering, o)
    , C(isotropicLinearElasticityTensor(E,v))
  {
    fs = apf::getShape(f);
  }
  void LinearElasticIntegrator::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    apf::NewArray<apf::Vector3> grads;
    apf::getShapeGrads(e,p,grads);
    // Linear strain-displacement matrix see Bathe pgs 555-556
    apf::DynamicMatrix B(6,nedofs);
    for(int ii = 0; ii < nenodes; ii++)
    {
      B(0,3*ii)   = grads[ii][0]; // N_(ii,1)
      B(0,3*ii+1) = B(0,3*ii+2) = 0.0;
      B(1,3*ii+1) = grads[ii][1]; // N_(ii,2)
      B(1,3*ii)   = B(1,3*ii+2) = 0.0;
      B(2,3*ii+2) = grads[ii][2]; // N_(ii,3)
      B(2,3*ii)   = B(2,3*ii+1) = 0.0;
      B(3,3*ii)   = grads[ii][1]; // N_(ii,2)
      B(3,3*ii+1) = grads[ii][0]; // N_(ii,1)
      B(3,3*ii+2) = 0.0;
      B(4,3*ii)   = 0.0;
      B(4,3*ii+1) = grads[ii][2]; // N_(ii,3)
      B(4,3*ii+2) = grads[ii][1]; // N_(ii,2)
      B(5,3*ii)   = grads[ii][2]; // N_(ii,3)
      B(5,3*ii+1) = 0.0;
      B(5,3*ii+2) = grads[ii][0]; // N_(ii,1)
    }
    apf::DynamicMatrix kt(nedofs,nedofs);
    apf::DynamicMatrix CB(6,nedofs);
    apf::multiply(C,B,CB);
    apf::DynamicMatrix BT(nedofs,6);
    apf::transpose(B,BT);
    apf::multiply(BT,CB,kt);
    // numerical integration
    kt *= w * dV;
    Ke += kt;
  }
}
