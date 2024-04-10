#include "apfElementalSystem.h"
#include <apfNumbering.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>
namespace amsi
{
  class apfElementalSystem : public ElementalSystem2
  {
  protected:
    apf::DynamicMatrix Ke;
    apf::DynamicVector Fe;
  public:
    apfElementalSystem(int nedofs)
      : ElementalSystem2(nedofs)
      , Ke(nedofs,nedofs)
      , Fe(nedofs)
    { zero(); }
    void zero()
    {
      Fe.zero();
      Ke.zero();
    }
    virtual double & fe(int idx)
    {
      return Fe(idx);
    }
    virtual double & ke(int rdx, int cdx)
    {
      return Ke(rdx,cdx);
    }
  };
  /*
   * creates a new elemental system
   * \warning destroyApfElementalSystem needs to be called after usage
   */
  ElementalSystem2 * buildApfElementalSystem(apf::Element * e, apf::Numbering * n)
  {
    int dof_per_nd = apf::countComponents(apf::getField(n));
    int nd_per_elm = apf::countNodes(e);
    int nedofs = dof_per_nd * nd_per_elm;
    ElementalSystem2 * es = new apfElementalSystem(nedofs);
    apf::NewArray<int> dofs;
    apf::getElementNumbers(n,apf::getMeshEntity(e),dofs);
    for(int ii = 0; ii < nedofs; ++ii)
      (*es).dofs(ii) = dofs[ii];
    return es;
  }
  void destroyApfElementalSystem(ElementalSystem2* es) {
    delete es;
  }
}
