#ifndef AMSI_BOUNDARY_CONDITIONS_H_
#define AMSI_BOUNDARY_CONDITIONS_H_
#include <amsiEnumOps.h>
#include <string>
namespace amsi
{
  #define BC_TYPES(OP) OP(dirichlet), OP(neumann), OP(num_bc_types)
  #define NEU_TYPES(OP) OP(custom), OP(traction), OP(pressure), OP(num_neu_types)
  enum BCType{BC_TYPES(MAKE_ENUM_OP)};
  enum NeuBCType{NEU_TYPES(MAKE_ENUM_OP)};
  int numBCComponents(int tp, int sbtp);
  int numDirichletComponents(int tp);
  int numNeumannComponents(int tp);
  struct BC
  {
    std::string fld;
    int tp;
    int sbtp;
  };
  class BCQuery
  {
  protected:
    BC * bc;
  public:
    BCQuery(BC * b)
      : bc(b)
    { }

  virtual bool isFixed(int ii = 0) = 0;
    virtual bool isConst(int ii = 0) = 0;
    virtual bool isTimeExpr(int ii = 0) = 0;
    virtual bool isSpaceExpr(int ii = 0) = 0;
    virtual double getValue(int ii = 0, ...) = 0;
  };
}
#endif
