#ifndef APF_MESH_UTIL_H_
#define APF_MESH_UTIL_H_
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
namespace amsi
{
  /**
   * Create an empty mesh associated with a 'null' model;
   * @return Pointer to an empty mesh.
   */
  apf::Mesh2 * makeNullMdlEmptyMesh();
  /**
   * Create a mesh containing a single entity related to a 'null' model.
   * @param t The entity type of the single mesh entity.
   * @param vs An array of correctly ordered coordinates defining the vertices of
   *           mesh entity to be created.
   * @return Pointer to a mesh with a single mesh entity.
   */
  apf::Mesh2 * makeSingleEntityMesh(apf::Mesh::Type t, const apf::Vector3 * vs);
  /**
   * Retrieve the first mesh entity given by iterating over the mesh. Most
   *  useful to retrieve the only mesh entity in a single entity mesh.
   */
  apf::MeshEntity * getFirstMeshEntity(apf::Mesh * msh, int dim);
  apf::Mesh2 * fiveTetCube();
}
#endif
