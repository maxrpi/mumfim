#ifndef AMSI_APF_BOUNDARY_CONDITIONS_H_
#define AMSI_APF_BOUNDARY_CONDITIONS_H_
#include "apfMesh.h"
#include <unordered_map>
namespace amsi
{
  typedef std::unordered_map<apf::MeshEntity*,std::vector<bool>> bc_set_map;
  class NeumannIntegrator;
  class BCQuery;
  class LAS;
  double getDirichletValue(BCQuery * qry,
                           apf::Mesh * msh,
                           apf::MeshEntity * ent,
                           int nd,
                           int cmp,
                           double t);
  template <typename I>
    int applyDirichletBC(apf::Numbering * nm,
                         const I& begin,
                         const I& end,
                         BCQuery * qry,
                         double t,
                         bc_set_map& already_set_numbering,
                         apf::Field* deltaField=NULL);
  template <typename I>
    void applyNeumannBC(LAS * las,
                        apf::Numbering * nm,
                        const I& bgn,
                        const I& nd,
                        NeumannIntegrator * i,
                        double t);
}
#include "apfBoundaryConditions_impl.h"
#endif
