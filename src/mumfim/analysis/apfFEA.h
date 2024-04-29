#ifndef AMSI_APFFEA_H_
#define AMSI_APFFEA_H_
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <iomanip>
#include "ElementalSystem.h"
#include "amsiBoundaryConditionQuery.h"
#include "amsiFEA.h"
#include "apfFieldOp.h"
namespace amsi {
  apf::Field* analyzeMeshQuality(apf::Mesh* mesh, apf::Field* disp_field);
  class PrintField : public amsi::FieldOp {
    private:
    apf::Field* f;
    apf::MeshEntity* me;
    std::ostream& os;
    int nc;
    public:
    PrintField(apf::Field* field, std::ostream& str)
        : f(field), me(), os(str), nc(0)
    {
      nc = apf::countComponents(f);
      os << std::setprecision(16);
    }
    virtual bool inEntity(apf::MeshEntity* e)
    {
      me = e;
      return true;
    }
    virtual void outEntity() {}
    virtual void atNode(int node)
    {
      double cmps[nc];
      apf::getComponents(f, me, node, &cmps[0]);
      for (int ii = 0; ii < nc; ii++)
        os << cmps[ii] << std::endl;
      os << std::endl;
    }
    void run() { apply(f); }
  };
}  // namespace amsi
#endif
