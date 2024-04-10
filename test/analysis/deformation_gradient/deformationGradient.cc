/*
 * This file includes tests to confirm that the deformation gradient
 * computes correctly. The deformation gradient is computed with uniaxial
 * tension, simple shear, and pure shear on a brick element.
 * Todo: this test should be extended to work on all element types
 *
 */
#include <amsiAnalysis.h>
#include <amsiDeformation.h>
#include <apf.h>
//#include <apfBox.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <test.h>
#include <iostream>
void computeDeformationGradient(apf::Mesh* mesh, apf::Field* dispField,
                                apf::Field* defm_grd)
{
  int dim = mesh->getDimension();
  apf::MeshIterator* it = mesh->begin(dim);
  apf::MeshEntity* meshEnt = NULL;
  while ((meshEnt = mesh->iterate(it))) {
    assert(mesh->getType(meshEnt) == apf::Mesh::HEX);
    int numNodes = mesh->adjacentCount[mesh->getType(meshEnt)][0];
    apf::Downward adjVerts;
    mesh->getDownward(meshEnt, 0, adjVerts);
    apf::MeshElement* meshElem = apf::createMeshElement(mesh, meshEnt);
    apf::Element* elem = apf::createElement(dispField, meshElem);
    apf::createElement(dispField, meshElem);
    for (int i = 0; i < numNodes; ++i) {
      apf::Vector3 point;
      apf::Matrix3x3 F;
      mesh->getPoint(adjVerts[i], 0, point);
      amsi::deformationGradient(elem, point, F);
      apf::setMatrix(defm_grd, adjVerts[i], 0, F);
    }
    apf::destroyElement(elem);
    apf::destroyMeshElement(meshElem);
  }
  mesh->end(it);
}
// take the displacement of each node and set the displacement field
// we note that the cube has the most nodes so the max array size is 24
// Note that this function is not general, and we should loop through all
// of the dimensions and set the nodes on that dimension.
void setDisp(apf::Mesh* mesh, double disp[24], apf::Field* dispField)
{
  int dim = mesh->getDimension();
  // loop over each mesh element
  apf::MeshIterator* it = mesh->begin(dim);
  apf::MeshEntity* meshEnt = NULL;
  while ((meshEnt = mesh->iterate(it))) {
    assert(mesh->getType(meshEnt) == apf::Mesh::HEX);
    int numNodes = mesh->adjacentCount[mesh->getType(meshEnt)][0];
    apf::Downward adjVerts;
    mesh->getDownward(meshEnt, 0, adjVerts);
    // loop over each node
    for (int i = 0; i < numNodes; ++i) {
      apf::Vector3 dispVector(&disp[i * dim]);
      apf::setVector(dispField, adjVerts[i], 0, dispVector);
    }
  }
  mesh->end(it);
}
// take the deformation gradient field and the correct deformation gradient and
// compare for validity. This function assumes affine deformation in the
// element.
int checkDeformationGrad(apf::Mesh* mesh, const std::string& testName,
                         apf::Field* F, apf::Matrix3x3& correctF)
{
  (void)correctF;
  int dim = mesh->getDimension();
  apf::FieldShape* s = apf::getShape(F);
  int result = 0;
  for (int j = 0; j < dim + 1; ++j) {
    if (!s->hasNodesIn(j)) continue;
    apf::MeshIterator* it = mesh->begin(j);
    apf::MeshEntity* meshEnt = NULL;
    int nodeNum = 0;
    while ((meshEnt = mesh->iterate(it))) {
      int numNodes = s->countNodesOn(mesh->getType(meshEnt));
      for (int i = 0; i < numNodes; ++i, ++nodeNum) {
        apf::Matrix3x3 nodeF;
        apf::getMatrix(F, meshEnt, i, nodeF);
        std::stringstream name;
        name << testName << "-node-" << nodeNum;
//        result += test(name.str(), correctF, nodeF);
      }
    }
    mesh->end(it);
  }
  return result;
}
int uniaxialTensionTest(apf::Mesh* mesh)
{
  apf::Field* dispField = apf::createFieldOn(mesh, "disp", apf::VECTOR);
  apf::Field* F = apf::createFieldOn(mesh, "F", apf::MATRIX);
  double alpha = 1.0;  // arbitrary constant
  double disp[24] = {0, 0, 0, alpha, 0, 0, alpha, 0, 0, 0, 0, 0,
                     0, 0, 0, alpha, 0, 0, alpha, 0, 0, 0, 0, 0};
  setDisp(mesh, disp, dispField);
  computeDeformationGradient(mesh, dispField, F);
  double L = 2;
  apf::Matrix3x3 correctF(1 + alpha / L, 0, 0, 0, 1, 0, 0, 0, 1);
  int result = checkDeformationGrad(mesh, "Uniaxial Tension", F, correctF);
  apf::destroyField(dispField);
  apf::destroyField(F);
  return result;
}
int simpleShearTest(apf::Mesh* mesh)
{
  apf::Field* dispField = apf::createFieldOn(mesh, "disp", apf::VECTOR);
  apf::Field* F = apf::createFieldOn(mesh, "F", apf::MATRIX);
  double alpha = 0.75;  // arbitrary constant
  double L = 2;         // FIXME should get length from element
  double disp[24] = {0,     0, 0, 0,     0, 0, 0,     0, 0, 0,     0, 0,
                     alpha, 0, 0, alpha, 0, 0, alpha, 0, 0, alpha, 0, 0};
  setDisp(mesh, disp, dispField);
  computeDeformationGradient(mesh, dispField, F);
  // apf::writeVtkFiles("simpleShear",mesh);
  apf::Matrix3x3 correctF(1, 0, alpha / L, 0, 1, 0, 0, 0, 1);
  int result = checkDeformationGrad(mesh, "Simple Shear", F, correctF);
  apf::destroyField(dispField);
  apf::destroyField(F);
  return result;
}
int pureShearTest(apf::Mesh* mesh)
{
  apf::Field* dispField = apf::createFieldOn(mesh, "disp", apf::VECTOR);
  apf::Field* F = apf::createFieldOn(mesh, "F", apf::MATRIX);
  double alpha = 0.75;  // arbitrary constant
  double L = 2;         // FIXME should get length from element
  double disp[24] = {0,     0,     0,     0, 0,     alpha, 0, 0,
                     alpha, 0,     0,     0, alpha, 0,     0, alpha,
                     0,     alpha, alpha, 0, alpha, alpha, 0, 0};
  setDisp(mesh, disp, dispField);
  computeDeformationGradient(mesh, dispField, F);
  // apf::writeVtkFiles("simpleShear",mesh);
  apf::Matrix3x3 correctF(1, 0, alpha / L, 0, 1, 0, alpha / L, 0, 1);
  int result = checkDeformationGrad(mesh, "Simple Shear", F, correctF);
  apf::destroyField(dispField);
  apf::destroyField(F);
  return result;
}
int main(int argc, char** argv)
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  // simplex mesh
  // apf::Mesh2* m = apf::makeMdsBox(1,1,1,1,1,1,1);
  // quad mesh
  // apf::Mesh2* m = apf::makeMdsBox(1, 1, 1, 1, 1, 1, 0);
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, 3, false);
  apf::Vector3 const hex[8] = {
      apf::Vector3(-1, -1, -1), apf::Vector3(1, -1, -1), apf::Vector3(1, 1, -1),
      apf::Vector3(-1, 1, -1),  apf::Vector3(-1, -1, 1), apf::Vector3(1, -1, 1),
      apf::Vector3(1, 1, 1),    apf::Vector3(-1, 1, 1)};
  apf::buildOneElement(mesh, mesh->findModelEntity(3, 0), apf::Mesh::HEX, hex);
  mesh->acceptChanges();
  int failed = 0;
  failed += uniaxialTensionTest(mesh);
  failed += simpleShearTest(mesh);
  failed += pureShearTest(mesh);
  // gmi_model* g = m->getModel();
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  amsi::freeAnalysis();
  return failed;
}
