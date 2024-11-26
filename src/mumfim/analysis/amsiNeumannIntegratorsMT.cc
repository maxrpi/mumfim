#include "amsiNeumannIntegratorsMT.h"
#include <array>

namespace amsi {
  NeumannIntegratorMT::NeumannIntegratorMT(LAS *l, apf::Field *f,
                                           const mt::IModelTrait *model_trait, int o,
                                           double t)
          : apf::Integrator(o)
          , las(l)
          , fe()
          , Ke()
          , fld(f)
          , me()
          , e()
          , nedofs()
          , nenodes()
          , nfcmps(apf::countComponents(fld))
          , tm(t)
          , mt_evaluator(model_trait)
  {
  }
  void NeumannIntegratorMT::setTime(double t) { tm = t; }
  int NeumannIntegratorMT::getnedofs() { return nedofs; }
  void NeumannIntegratorMT::inElement(apf::MeshElement *m)
  {
    me = m;
    e = apf::createElement(fld, me);
    nenodes = apf::countNodes(e);
    int nnedofs = nenodes * nfcmps;
    if (nnedofs != nedofs)
    {
      fe.setSize(nnedofs);
      Ke.setSize(nnedofs, nnedofs);
    }
    nedofs = nnedofs;
    fe.zero();
    Ke.zero();
  }
  void NeumannIntegratorMT::outElement() { apf::destroyElement(e); }
  SurfaceTractionMT::SurfaceTractionMT(LAS *l, apf::Field *f,
                                       const mt::IModelTrait *mt, int o, double t)
      : NeumannIntegratorMT(l, f, mt, o, t)
  {
  }
  void SurfaceTractionMT::atPoint(const apf::Vector3 &p, double w, double dV)
  {
    // mt->visit(visitor);
    // loads the value for t, p(xyz point) into vls
    std::vector<double> vals;
    apf::Vector3 xyz;
    apf::mapLocalToGlobal(me, p, xyz);
    mt_evaluator(tm, xyz[0], xyz[1], xyz[2], vals);
    if (vals.size() != static_cast<size_t>(nfcmps)) {
      std::cerr << "The Surface Traction model trait must have the same number of "
                   "components as the field\n";
      exit(1);
    }
    apf::NewArray<double> N;
    apf::getShapeValues(e, p, N);
    double wxdV = w * dV;
    for (int nd = 0; nd < nenodes; nd++) {
      for (int cmp = 0; cmp < nfcmps; cmp++) {
        fe(nd * nfcmps + cmp) = N[nd] * vals[cmp] * wxdV;
      }
    }
  }

  PressureMT::PressureMT(LAS *l, apf::Field *f, const mt::IModelTrait *mt, int o,
                         double t)
      : NeumannIntegratorMT(l, f, mt, o, t), msh(apf::getMesh(f))
  {
  }

  void PressureMT::inElement(apf::MeshElement *m)
  {
    NeumannIntegratorMT::inElement(m);
    ent = apf::getMeshEntity(m);
  }

  void PressureMT::atPoint(const apf::Vector3 &p, double w, double dV)
  {
    std::vector<double> vals;
    apf::Vector3 xyz;
    apf::mapLocalToGlobal(me, p, xyz);
    mt_evaluator(tm, xyz[0], xyz[1], xyz[2], vals);
    if (vals.size() != 1) {
      std::cerr << "The Pressure model trait is expected to have 1 component.\n";
      exit(1);
    }
    apf::NewArray<double> N;
    std::array<double,3> pressure;
    apf::getShapeValues(e, p, N);
    apf::Vector3 nrml;
    faceNormal(msh, ent, nrml);
    pressure[0] = nrml.x()*vals[0];
    pressure[1] = nrml.y()*vals[0];
    pressure[2] = nrml.z()*vals[0];
    double wxdV = w * dV;
    for (int nd = 0; nd < nenodes; nd++) {
      for (int cmp= 0; cmp < nfcmps; cmp++) {
        fe(nd * nfcmps + cmp) = N[nd] * vals[cmp] * wxdV;
      }
    }
  }

  RobinMT::RobinMT(LAS *l, apf::Field *f, const mt::IModelTrait *mt, int o,
                         double t)
      : NeumannIntegratorMT(l, f, mt, o, t), msh(apf::getMesh(f))
  {
    if (nfcmps != 1) {
      std::cerr << "The Robin BC is only defined for scalar fields.\n";
      exit(1);
    }
  }

  void RobinMT::inElement(apf::MeshElement *m)
  {
    NeumannIntegratorMT::inElement(m);
    ent = apf::getMeshEntity(m);
  }

  void RobinMT::atPoint(const apf::Vector3 &p, double w, double dV)
  {
    std::vector<double> vals;
    apf::Vector3 xyz;
    apf::mapLocalToGlobal(me, p, xyz);
    mt_evaluator(tm, xyz[0], xyz[1], xyz[2], vals);
    if (vals.size() != 2) {
      std::cerr << "The Robin model trait is expected to have 2 parameter values (h, T_inf).\n";
      exit(1);
    }

    // This is likely for convective heat transfer, but could be for mass transfer, etc.
    // For thermal problems transfer_coefficient = heat transfer coefficient,
    // and field_at_inf = temperature across the boundary
    double transfer_coefficient = vals[0];
    double field_at_inf = vals[1];
    apf::NewArray<double> N;
    apf::getShapeValues(e, p, N);
    for (int i = 0; i < nenodes; i++){
      fe(i) = N[i] * transfer_coefficient * field_at_inf * w * dV;
      for (int j = 0; j < nenodes; j++)
      {
        Ke(i,j) = N[i] * N[j] * transfer_coefficient * w * dV;
      }
    }
  }


  std::unique_ptr<NeumannIntegratorMT> createNeumannIntegrator(
      LAS *las, apf::Field *fld, const mt::IModelTrait *mt, int o, double t,
      NeumannBCType tp)
  {
    switch(tp)
    {
      case NeumannBCType::traction:
        return std::make_unique<SurfaceTractionMT>(las,fld,mt, o,t);
      case NeumannBCType::pressure:
        return std::make_unique<PressureMT>(las,fld,mt,o,t);
      case NeumannBCType::robin:
        return std::make_unique<RobinMT>(las,fld,mt,o,t);
    }
    return {nullptr};
  }
}
