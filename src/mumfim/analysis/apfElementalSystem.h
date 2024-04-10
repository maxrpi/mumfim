#ifndef APF_ELEMENTAL_SYSTEM_H_
#define APF_ELEMENTAL_SYSTEM_H_
#include "amsiElementalSystem.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfDynamicMatrix.h>
namespace amsi
{
  ElementalSystem2 * buildApfElementalSystem(apf::Element * e, apf::Numbering * n);
  void destroyApfElementalSystem(ElementalSystem2* es);
}
#endif
