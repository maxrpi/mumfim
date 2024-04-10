#include "amsiBoundaryConditionQuery.h"
#include <cassert>
#include "amsiFields.h"
namespace amsi
{
  const char * const BCTypes[] = {BC_TYPES(MAKE_STRING_OP)};
  const char * const NeuTypes[] = {NEU_TYPES(MAKE_STRING_OP)};
  char const * getBCTypeString(int tp)
  {
    assert(tp < BCType::num_bc_types);
    return BCTypes[tp];
  }
  char const * getBCSubtypeString(int tp, int sbtp)
  {
    switch(tp)
    {
    case BCType::dirichlet:
      return fieldUnitString(static_cast<FieldUnit>(sbtp));
    case BCType::neumann:
      return getNeumannTypeString(static_cast<FieldType>(sbtp));
    default:
      return NULL;
    }
  }
  char const * getNeumannTypeString(int tp)
  {
    assert(tp < NeuBCType::num_neu_types && tp >= 0);
    return NeuTypes[tp];
  }
  int numBCComponents(int tp, int sbtp)
  {
    if(tp == BCType::dirichlet)
      return numDirichletComponents(sbtp);
    else if(tp == BCType::neumann)
      return numNeumannComponents(sbtp);
    else
      return 0;
  }
  int numDirichletComponents(int tp)
  {
    switch(tp)
    {
    case FieldUnit::unitless:
    case FieldUnit::displacement:
      return 3;
    default:
      return 0;
    }
  }
  int numNeumannComponents(int tp)
  {
    switch(tp)
    {
    case NeuBCType::traction:
      return 3;
    case NeuBCType::pressure:
      return 1;
    default:
      return 0;
    }
  }
}
