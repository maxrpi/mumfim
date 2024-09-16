#ifndef MUMFIM_FOURIER_INTEGRATOR_H_
#define MUMFIM_FOURIER_INTEGRATOR_H_
#include "ElementalSystem.h"

namespace mumfim{

  class EffectiveKappaIntegrator : public amsi::ElementalSystem
  {
  apf::DynamicMatrix D;  // elemental kappa
  apf::MeshElement *mesh_element;
  std::map<std::pair<int, int>, double *> &k_map;
  apf::Matrix3x3 * getIsotropicKappa(double kii);
  apf::Matrix3x3 * getFullKappa(double xx, double yy, double zz,
                                  double xy, double yz, double zx);
  public:
    EffectiveKappaIntegrator(apf::Field *, apf::Numbering* numbering, apf::Matrix3x3 * K_r, std::map<std::pair<int, int>, double *> &k_map_);
    void atPoint(apf::Vector3 const &p,
      double w, double dV);

  protected:
  };
}
#endif