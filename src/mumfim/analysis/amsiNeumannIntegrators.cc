#include "amsiNeumannIntegrators.h"
#include "amsiBoundaryConditionQuery.h"
namespace amsi
{
  void NeumannIntegrator::updateBCQueryValues(apf::Vector3 const & p)
  {
    // nfcmps is the number of components on the field
    // this assumes that the the query has the same number of
    // components as the field...should probably be tested
    if(qry->isConst())
      for(int ii = 0; ii < nfcmps; ii++)
        vls[ii] = qry->getValue(ii);
    else if(qry->isSpaceExpr())
    {
      apf::Vector3 xyz;
      apf::mapLocalToGlobal(me,p,xyz);
      if(qry->isTimeExpr())
        for(int ii = 0; ii < nfcmps; ii++)
          vls[ii] = qry->getValue(ii,tm,xyz[0],xyz[1],xyz[2]);
      else
        for(int ii = 0; ii < nfcmps; ii++)
          vls[ii] = qry->getValue(ii,xyz[0],xyz[1],xyz[2]);
    }
    else if(qry->isTimeExpr())
      for(int ii = 0; ii < nfcmps; ii++)
        vls[ii] = qry->getValue(ii,tm);
  }
  NeumannIntegrator * buildNeumannIntegrator(LAS * las, apf::Field * fld, int o, BCQuery * qry, int tp, double t)
  {
    switch(tp)
    {
    case NeuBCType::traction:
      return new SurfaceTraction(las,fld,o,qry,t);
    case NeuBCType::pressure:
      return new Pressure(las,fld,o,qry,t);
    default:
      return NULL;
    }
  }
  void deleteNeumannIntegrator(NeumannIntegrator * i)
  {
    delete i;
  }
}
