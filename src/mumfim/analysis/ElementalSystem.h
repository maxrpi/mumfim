#ifndef ELEMENTALSYSTEM_H_
#define ELEMENTALSYSTEM_H_
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfNumbering.h>
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
    apf::Numbering* numbering_;
    int nedofs;
    int nenodes;
    int num_field_components;
    std::vector<double> field_values_;
    std::vector<int> field_numbers_;
  public:
    ElementalSystem(apf::Field * f, apf::Numbering* numbering, int o);
    virtual void inElement(apf::MeshElement *);
    virtual void outElement();
    virtual void parallelReduce() {};
    virtual bool includesBodyForces() {return false;}
    virtual int numElementalDOFs() {return nedofs;}
    virtual apf::DynamicMatrix& getKe() {return Ke;}
    virtual apf::DynamicVector& getfe() {return fe;}
    virtual void zeroKe(){Ke.zero();}
    virtual void zerofe(){fe.zero();}

    // TODO have getFieldValues and getFieldNumbers return mdspan
    [[nodiscard]] const std::vector<double>& getFieldValues() { return field_values_; }
    [[ nodiscard ]] const std::vector<int> & getFieldNumbers() { return field_numbers_; }
    virtual apf::Field * getField() {return f;}
  };
}
#endif
