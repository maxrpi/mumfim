#ifndef AMSI_NEUMANN_INTEGRATORS_H_
#define AMSI_NEUMANN_INTEGRATORS_H_
#include "amsiLAS.h"
#include "apfFunctions.h"
#include <apf.h>
#include <apfDynamicVector.h>
#include <cassert>
namespace amsi
{
  class BCQuery;
  class NeumannIntegrator : public apf::Integrator
  {
  protected:
    LAS * las;
    apf::DynamicVector fe;
    apf::Field * fld;
    apf::MeshElement * me;
    apf::Element * e;
    int nedofs;
    int nenodes;
    int nfcmps;
    double tm;
    std::vector<int> dofs;
    std::vector<double> vls;
    BCQuery * qry;
  public:
  NeumannIntegrator(LAS * l, apf::Field * f, int o, BCQuery * q, double t = 0.0)
      : apf::Integrator(o)
      , las(l)
      , fe()
      , fld(f)
      , me()
      , e()
      , nedofs()
      , nenodes()
      , nfcmps(apf::countComponents(fld))
      , tm(t)
      , vls(nfcmps)
      , qry(q)
    { }
    void setTime(double t) { tm = t; }
    int getnedofs() { return nedofs; }
    apf::DynamicVector & getFe() { return fe; }
    virtual void inElement(apf::MeshElement * m)
    {
      me = m;
      e = apf::createElement(fld,me);
      nenodes = apf::countNodes(e);
      //assert(nfcmps == qry->numComps());
      int nnedofs = nenodes * nfcmps;
      if(nnedofs != nedofs)
        fe.setSize(nnedofs);
      nedofs = nnedofs;
      fe.zero();
    }
    virtual void outElement()
    {
      apf::destroyElement(e);
    }
  protected:
    void updateBCQueryValues(apf::Vector3 const & p);
  };
  class SurfaceTraction : public NeumannIntegrator
  {
  public:
    SurfaceTraction(LAS * l, apf::Field * f, int o, BCQuery * q, double t = 0.0)
      : NeumannIntegrator(l,f,o,q,t)
    { }
    void atPoint(apf::Vector3 const & p, double w, double dV)
    {
      // loads the value for t, p(xyz point) into vls
      updateBCQueryValues(p);
      apf::NewArray<double> N;
      apf::getShapeValues(e,p,N);
      double wxdV = w * dV;
      for(int ii = 0; ii < nenodes; ii++)
        for(int jj = 0; jj < nfcmps; jj++)
          fe(ii*nfcmps + jj) = N[ii] * vls[jj] * wxdV;
    }
  };
  class Pressure : public NeumannIntegrator
  {
  private:
    apf::Mesh * msh;
    apf::MeshEntity * ent;
  public:
    Pressure(LAS * l, apf::Field * f, int o, BCQuery * q, double t)
      : NeumannIntegrator(l,f,o,q,t)
      , msh(apf::getMesh(f))
    { }
    void inElement(apf::MeshElement * m)
    {
      NeumannIntegrator::inElement(m);
      ent = apf::getMeshEntity(m);
    }
    void atPoint(apf::Vector3 const & p, double w, double dV)
    {
      updateBCQueryValues(p);
      apf::NewArray<double> N;
      apf::getShapeValues(e,p,N);
      apf::Vector3 nrml;
      faceNormal(msh,ent,nrml);
      vls[0] *= nrml.x();
      vls[1] *= nrml.y();
      vls[2] *= nrml.z();
      double wxdV = w * dV;
      for(int ii = 0; ii < nenodes; ii++)
        for(int jj = 0; jj < nfcmps; jj++)
          fe(ii*nfcmps + jj) = N[ii] * vls[jj] * wxdV;
    }
  };
  NeumannIntegrator * buildNeumannIntegrator(LAS * las,
                                             apf::Field * fld,
                                             int o,
                                             BCQuery * qry,
                                             int tp,
                                             double t = 0.0);
  void deleteNeumannIntegrator(NeumannIntegrator * i);
}
#endif
