#ifndef AMSI_ELEMENTAL_SYSTEM_H_
#define AMSI_ELEMENTAL_SYSTEM_H_
namespace amsi
{
  class ElementalSystem2
  {
  protected:
    int nedofs;
    int * edofs;
  public:
    ElementalSystem2(int n)
      : nedofs(n)
      , edofs(new int [nedofs])
    { }
    virtual ~ElementalSystem2()
    {
      delete [] edofs;
    }
    virtual void zero() = 0;
    virtual double & fe(int idx) = 0;
    virtual double & ke(int rdx, int cdx) = 0;
    int & dofs(int ddx)
    {
      return edofs[ddx];
    }
    int nedof()
    {
      return nedofs;
    }
  };
}
#endif
