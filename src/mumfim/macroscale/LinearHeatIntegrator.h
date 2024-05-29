#ifndef MUMFIM_FOURIER_INTEGRATOR_H_
#define MUMFIM_FOURIER_INTEGRATOR_H_
#include "ElementalSystem.h"

namespace mumfim{

  class LinearHeatIntegrator : public amsi::ElementalSystem
  {
  apf::Field * T_;
  apf::Field * K_;
  apf::Field * f_;
  apf::DynamicMatrix D;  // elemental kappa
  apf::Matrix3x3 * getIsotropicKappa(double kii);
  apf::Matrix3x3 * getFullKappa(double xx, double yy, double zz,
                                  double xy, double yz, double zx);
  public:
    LinearHeatIntegrator(apf::Field *T, apf::Numbering* numbering, apf::Field * K);
    LinearHeatIntegrator(apf::Field *T, apf::Numbering* numbering, apf::Matrix3x3 * K_r); // Piecewise constant Kappa
    //void inElement(apf::MeshElement *me);
    bool includesBodyForces() final { return true; }
    void atPoint(apf::Vector3 const &p,
      double w, double dV);

  protected:
  };
}
#endif