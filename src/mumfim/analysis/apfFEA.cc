#include "apfFEA.h"
#include <apfNumbering.h>
#include <apfShape.h>
#include <maShape.h>
#include <cassert>
#include "amsiBoundaryConditions.h"
#include "apfFunctions.h"
namespace amsi {
  // assumes only linear tets
  apf::Field* analyzeMeshQuality(apf::Mesh* mesh, apf::Field* disp_field)
  {
    int analysis_dim = mesh->getDimension();
    apf::Field* f = createStepField(mesh, "quality", apf::SCALAR);
    apf::MeshEntity* me = NULL;
    for (apf::MeshIterator* it = mesh->begin(analysis_dim);
         (me = mesh->iterate(it));) {
      apf::Vector3 verts[4];
      apf::Adjacent v;
      mesh->getAdjacent(me, 0, v);
      for (int ii = 0; ii < 4; ii++) {
        apf::Vector3 disp;
        apf::getVector(disp_field, v[ii], 0, disp);
        mesh->getPoint(v[ii], 0, verts[ii]);
        verts[ii] += disp;
      }
      apf::setScalar(f, me, 0, ma::measureLinearTetQuality(verts));
    }
    return f;
  }
}  // namespace amsi
