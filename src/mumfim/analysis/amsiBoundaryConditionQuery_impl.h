#ifndef AMSI_BOUNDARY_CONDITIONS_IMPL_H_
#define AMSI_BOUNDARY_CONDITIONS_IMPL_H_
#include "amsiFields.h"
namespace amsi
{
  // probably need a richer set of tools to manage the supported field types and their associated boundary conditions (ie a metamodel)
  template <typename O>
    void getApplicableBCTypesForField(int fld_tp, int bc_tp, O out)
  {
    if(fld_tp == FieldUnit::displacement && bc_tp == BCType::dirichlet)
      *out++ = displacement;
    else if(fld_tp == FieldUnit::displacement && bc_tp == BCType::neumann)
    {
      *out++ = NeuBCType::traction;
      *out++ = NeuBCType::pressure;
    }
  }
}
#endif
