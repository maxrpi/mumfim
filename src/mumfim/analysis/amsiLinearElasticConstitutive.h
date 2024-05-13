#ifndef AMSI_LINEAR_ELASTIC_CONSTITUTIVE_H_
#define AMSI_LINEAR_ELASTIC_CONSTITUTIVE_H_
#include "ElementalSystem.h"
namespace amsi
{
  apf::DynamicMatrix isotropicLinearElasticityTensor(double E, double v);
  class LinearElasticIntegrator : public ElementalSystem
  {
  public:
    LinearElasticIntegrator(apf::Field * field,
                            apf::Numbering* numbering,
                            int o,
                            double E,
                            double v);
    void atPoint(apf::Vector3 const &p,
                 double w,
                 double dV);
  private:
    apf::FieldShape * fs;
    apf::DynamicMatrix C;
  };
}
#endif
