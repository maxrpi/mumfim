#ifndef MUMFIM_FOURIER_INTEGRATOR_H_
#define MUMFIM_FOURIER_INTEGRATOR_H_
#include "ElementalSystem.h"

namespace mumfim{

  // Only integrates the stiffness matrix for linear heat conduction. Nothing from the RHS.
  class ThermalStiffnessIntegrator : public amsi::ElementalSystem
  {
  apf::Field * f_;
  apf::DynamicMatrix D;  // elemental kappa
  apf::Matrix3x3 * getIsotropicKappa(double kii);
  apf::Matrix3x3 * getFullKappa(double xx, double yy, double zz,
                                  double xy, double yz, double zx);
  public:
    ThermalStiffnessIntegrator(apf::Field *T, apf::Numbering* numbering, apf::Field * K);
    ThermalStiffnessIntegrator(apf::Field *T, apf::Numbering* numbering, apf::Matrix3x3 * K_r); // Piecewise constant Kappa
    //void inElement(apf::MeshElement *me);
    void atPoint(apf::Vector3 const &p,
      double w, double dV);

  protected:
  };
}
#endif