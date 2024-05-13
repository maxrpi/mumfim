#include "ElementalSystem.h"
#include <apfNumbering.h>
#include <apfField.h>

namespace amsi
{
  // this function works for both scalar and vector fields
  static void GetNodalFieldValuesAndNumbers(apf::Field* field,
                                            apf::Numbering* numbering,
                                            apf::MeshEntity* mesh_entity, std::vector<int>& dof_numbers, std::vector<double> & values) {
    auto nnodes = field->countNodesOn(mesh_entity);
    const auto ncomponents = apf::countComponents(field);
    values.resize(ncomponents*nnodes);
    dof_numbers.resize(ncomponents*nnodes);
    for(int node=0; node<nnodes; ++node) {
      apf::getComponents(field, mesh_entity, node, &values[node*ncomponents]);
      for(int component=0; component<ncomponents; ++component) {
        dof_numbers[node*ncomponents+component] = apf::getNumber(numbering,mesh_entity,node, component);
      }
    }
  }
  ElementalSystem::ElementalSystem(apf::Field * field, apf::Numbering* numbering, int o)
      : apf::Integrator(o)
      , Ke()
      , fe()
      , me(NULL)
      , e(NULL)
      , f(field)
      , numbering_(numbering)
      , nedofs(0)
      , nenodes(0)
      , num_field_components(0)
  {
    num_field_components = apf::countComponents(f);
  }

  // if the new mesh elemenent is on a mesh entity with the same type as the
  // previous mesh element was, we don't need to reallocate the memory for the
  // element matrix and vector (just zero it), otherwise we need to resize it
  void ElementalSystem::inElement(apf::MeshElement * ME)
  {
    me = ME;
    e = apf::createElement(f, me);
    nenodes = apf::countNodes(e);
    int new_nedofs = nenodes * num_field_components;
    bool reallocate = nedofs != new_nedofs;
    nedofs = new_nedofs;
    GetNodalFieldValuesAndNumbers(f, numbering_, getMeshEntity(ME),
                                  field_numbers_, field_values_);
        // TODO: remove this check. login in DynamicArray will account for this
    // Note: the reallocation routine in DynamicArray will cause a lot of reallocations
    if (reallocate)
    {
      Ke.setSize(nedofs, nedofs);
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
}  // namespace amsi
