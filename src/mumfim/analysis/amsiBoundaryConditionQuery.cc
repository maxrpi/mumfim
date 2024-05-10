#include "amsiBoundaryConditionQuery.h"
#include <cassert>
#include "amsiFields.h"
namespace amsi
{
  const char * const BCTypes[] = {BC_TYPES(MAKE_STRING_OP)};
  const char * const NeuTypes[] = {NEU_TYPES(MAKE_STRING_OP)};

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
