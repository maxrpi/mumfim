#include "ElementalSystem.h"
namespace amsi
{
  ElementalSystem::ElementalSystem(apf::Field * field, int o)
    : apf::Integrator(o)
    , Ke()
    , fe()
    , me(NULL)
    , e(NULL)
    , f(field)
    , nedofs(0)
    , nenodes(0)
    , num_field_components(0)
  {
    num_field_components = apf::countComponents(f);
  }
  // if the new mesh elemenent is on a mesh entity with the same type as the previous mesh element was, we don't need to reallocate the memory for the element matrix and vector (just zero it), otherwise we need to resize it
  void ElementalSystem::inElement(apf::MeshElement * ME)
  {
    e = apf::createElement(f,me);
    nenodes = apf::countNodes(e);
    int new_nedofs = nenodes * num_field_components;
    bool reallocate = nedofs != new_nedofs;
    nedofs = new_nedofs;
    if(reallocate)
    {
      Ke.setSize(nedofs,nedofs);
      fe.setSize(nedofs);
    }
    Ke.zero();
    fe.zero();
  }
  void ElementalSystem::outElement()
  {
    apf::destroyElement(e);
    e = nullptr;
  }
}
