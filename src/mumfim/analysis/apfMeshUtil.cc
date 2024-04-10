#include "apfMeshUtil.h"
#include <gmi.h>
#include <gmi_null.h>
#include <apfMDS.h>
namespace amsi
{
  apf::Mesh2 * makeNullMdlEmptyMesh()
  {
    gmi_register_null();
    gmi_model * mdl = gmi_load(".null");
    return apf::makeEmptyMdsMesh(mdl,3,false);
  };
  apf::Mesh2 * makeSingleEntityMesh(apf::Mesh::Type t, const apf::Vector3 * vs)
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    apf::buildOneElement(msh,NULL,t,vs);
    // only required in parallel
    msh->acceptChanges();
    if(t != apf::Mesh::EDGE)
      apf::reorderMdsMesh(msh);
    apf::deriveMdsModel(msh);
    return msh;
  }
  apf::MeshEntity * getFirstMeshEntity(apf::Mesh * msh, int d)
  {
    apf::MeshEntity * ent = NULL;
    apf::MeshIterator * itr = msh->begin(d);
    ent = msh->iterate(itr);
    msh->end(itr);
    return ent;
  }
  const apf::Vector3 five_tet_vrts[8] =
  {
    apf::Vector3(-1.0,-1.0, 1.0),
    apf::Vector3( 1.0,-1.0, 1.0),
    apf::Vector3( 1.0,-1.0,-1.0),
    apf::Vector3(-1.0,-1.0,-1.0),
    apf::Vector3(-1.0, 1.0, 1.0),
    apf::Vector3( 1.0, 1.0, 1.0),
    apf::Vector3( 1.0, 1.0,-1.0),
    apf::Vector3(-1.0, 1.0,-1.0)
  };
  const int five_tet_order[5][4] =
  {
    {0,1,3,4},
    {1,3,2,6},
    {1,3,4,6},
    {4,6,5,1},
    {4,6,7,2}
  };
  apf::Mesh2 * fiveTetCube()
  {
    apf::Mesh2 * msh = makeNullMdlEmptyMesh();
    for(int tt = 0; tt < 5; ++tt)
    {
      apf::MeshEntity * vrts[4];
      for(int vv = 0; vv < 4; ++vv)
        vrts[vv] = msh->createVertex(NULL,five_tet_vrts[five_tet_order[tt][vv]],apf::Vector3(0.0,0.0,0.0));
      apf::buildElement(msh,NULL,apf::Mesh::TET,&vrts[0]);
    }
    return msh;
  }
}
