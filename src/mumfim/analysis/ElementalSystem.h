#ifndef ELEMENTALSYSTEM_H_
#define ELEMENTALSYSTEM_H_
#include <apf.h>
#include <apfDynamicMatrix.h>
namespace amsi
{
  class ElementalSystem : public apf::Integrator
  {
  protected:
    apf::DynamicMatrix Ke;
    apf::DynamicVector fe;
    apf::MeshElement * me;
    apf::Element * e;
    apf::Field * f;
    int nedofs;
    int nenodes;
    int num_field_components;
    virtual void _inElement(apf::MeshElement *) {}
  public:
    ElementalSystem(apf::Field * f, int o);
    virtual void inElement(apf::MeshElement *);
    virtual void outElement();
    virtual void parallelReduce() {};
    virtual bool includesBodyForces() {return false;}
    virtual int numElementalDOFs() {return nedofs;}
    // This shouldn't be used because it causes segfaults when
    // used improperly (which is almost always!
    // virtual apf::Element * getElement() {return e;}
    virtual apf::DynamicMatrix& getKe() {return Ke;}
    virtual apf::DynamicVector& getfe() {return fe;}
    virtual apf::Field * getField() {return f;}
  };
}
#endif
