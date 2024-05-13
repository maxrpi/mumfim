#include "apfMeasure.h"
#include <amsiMPI.h>
#include <apfShape.h>
#include "ElementalSystem.h"
namespace amsi {
  double distance(const apf::Vector3& pt_a, const apf::Vector3& pt_b)
  {
    double len = 0.0;
    for (int ii = 0; ii < 3; ii++)
      len += (pt_a[ii] - pt_b[ii]) * (pt_a[ii] - pt_b[ii]);
    return std::sqrt(len);
  }
  double triangleArea(const apf::Vector3& pt_a, const apf::Vector3& pt_b,
                      const apf::Vector3& pt_c)
  {
    // Find area of triangle from 3 points in 3D space using Heron's formula
    double len_a = distance(pt_b, pt_c);
    double len_b = distance(pt_a, pt_c);
    double len_c = distance(pt_a, pt_b);
    double s = 0.5 * (len_a + len_b + len_c);
    return std::sqrt(s * (s - len_a) * (s - len_b) * (s - len_c));
  }
  void faceNormal(const apf::Vector3& pt_a, const apf::Vector3& pt_b,
                  const apf::Vector3& pt_c, apf::Vector3& n)
  {
    apf::Plane p = apf::Plane::fromPoints(pt_a, pt_b, pt_c);
    n = p.normal;
  }
  int side(apf::ModelEntity* mdl_ent, apf::Mesh* msh, apf::MeshEntity* fc)
  {
    int nrml = 0;
    int tg = msh->getModelTag(mdl_ent);
    apf::MeshEntity* rgn0 = msh->getUpward(fc, 0);
    apf::MeshEntity* rgn1 = msh->getUpward(fc, 1);
    int tg0 = rgn0 == nullptr ? -2 : msh->getModelTag(msh->toModel(rgn0));
    int tg1 = rgn1 == nullptr ? -2 : msh->getModelTag(msh->toModel(rgn1));
    if (tg0 == tg1)
      nrml = 0;
    else if (tg0 == tg)
      nrml = 1;
    else if (tg1 == tg)
      nrml = -1;
    else  // how can this ever happen?!?!
    {
      if (tg0 == -2 && tg1 != tg)
        nrml = 1;
      else if (tg1 == -2 && tg0 != tg)
        nrml = -1;
    }
    return nrml;
  }
  class MeasureDisplaced : public amsi::ElementalSystem {
    private:
    int dim;
    int ent_dim;
    apf::FieldShape* fs;
    apf::EntityShape* es;
    apf::Element* mesh_coord_elem;
    double vol;

    public:
    MeasureDisplaced(apf::Field* field, apf::Numbering* numbering, int o)
        : ElementalSystem(field, numbering, o)
        , dim(apf::getMesh(field)->getDimension())
        , ent_dim(-1)
        , fs(NULL)
        , es(NULL)
        , mesh_coord_elem(NULL)
        , vol(0)
    {
    }
    void inElement(apf::MeshElement* me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      ent_dim = apf::getDimension(me);
      mesh_coord_elem =
          apf::createElement(apf::getMesh(f)->getCoordinateField(), me);
    }
    void atPoint(apf::Vector3 const&, double w, double)
    {
      int& nen = nenodes;  // = 4 (tets)
      // 1. Get coordinates on underlying mesh
      apf::NewArray<apf::Vector3> mesh_xyz;
      apf::getVectorNodes(mesh_coord_elem, mesh_xyz);
      // 2. Get coordinates from apf_primary_field (passed in), which contains
      // the accumulated displacement
      apf::NewArray<apf::Vector3> primary_field_xyz;
      apf::getVectorNodes(e, primary_field_xyz);
      // 3. Calculate current coordinates
      apf::DynamicMatrix xyz(nen, dim);
      xyz.zero();
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          xyz(ii, jj) = mesh_xyz[ii][jj] + primary_field_xyz[ii][jj];
      // For Updated Lagrangian, the Jacobian of the updated coordinates are
      // used Note: that entires of Jacobian is hard coded for Linear tetrahedra
      // elements. TO DO: Generalize Jacobian for current configuration.
      apf::Matrix<3, 3> J;
      J[0][0] = xyz(1, 0) - xyz(0, 0);  // x2-x1
      J[0][1] = xyz(2, 0) - xyz(0, 0);  // x3-x1
      J[0][2] = xyz(3, 0) - xyz(0, 0);  // x4-x1
      J[1][0] = xyz(1, 1) - xyz(0, 1);  // y2-y1
      J[1][1] = xyz(2, 1) - xyz(0, 1);  // y3-y1
      J[1][2] = xyz(3, 1) - xyz(0, 1);  // y4-y1
      J[2][0] = xyz(1, 2) - xyz(0, 2);  // z2-z1
      J[2][1] = xyz(2, 2) - xyz(0, 2);  // z3-z1
      J[2][2] = xyz(3, 2) - xyz(0, 2);  // z4-z1
      double detJ = getDeterminant(J);
      vol += w * detJ;
    }
    void outElement()
    {
      apf::destroyElement(mesh_coord_elem);
      ElementalSystem::outElement();
    }
    double getVol() { return vol; }
  };
  class MeasureDisplacedFromSurf : public amsi::ElementalSystem {
    public:
    MeasureDisplacedFromSurf(apf::Field* field, apf::Numbering* numbering, int o, int normal_dir)
        : ElementalSystem(field, numbering, o), dim(0), vol(0), norm_dir(normal_dir)
    {
    }
    void inElement(apf::MeshElement* me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      // dim = apf::getDimension(me);
      /** We want to consider 3D space, however MeshElement me is for 2D entity.
       * Therefore we need to hardcode dim for now.. */
      dim = 3;
    }
    void atPoint(apf::Vector3 const&, double, double)
    {
      int& nen = nenodes;  // = 3 (triangle)
      // 1. Get coordinates on underlying mesh
      apf::Mesh* mesh = apf::getMesh(f);
      apf::Field* apf_coord_field = mesh->getCoordinateField();
      apf::Element* mesh_coord_elem = apf::createElement(apf_coord_field, me);
      apf::NewArray<apf::Vector3> mesh_xyz;
      apf::getVectorNodes(mesh_coord_elem, mesh_xyz);
      // 2. Get coordinates from apf_primary_field (passed in), which contains
      // the accumulated displacement
      apf::NewArray<apf::Vector3> primary_field_xyz;
      apf::getVectorNodes(e, primary_field_xyz);
      // 3. Calculate current coordinates
      apf::DynamicMatrix xyz(nen, dim);
      xyz.zero();
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          xyz(ii, jj) = mesh_xyz[ii][jj] + primary_field_xyz[ii][jj];
      // Calculate volumes from triangles on surface mesh
      apf::Vector3 normal;
      apf::Vector3 pt0, pt1, pt2;
      for (int jj = 0; jj < dim; jj++) {
        pt0[jj] = xyz(0, jj);
        pt1[jj] = xyz(1, jj);
        pt2[jj] = xyz(2, jj);
      }
      double area = triangleArea(pt0, pt1, pt2);
      faceNormal(pt0, pt1, pt2, normal);
      vol = 0.0;
      for (int jj = 0; jj < dim; jj++)
        vol += norm_dir * normal[jj] * (pt0[jj] + pt1[jj] + pt2[jj]);
      vol *= area / 3.0;
      vol *= 1.0 / 3.0;
    }
    double getVol() { return vol; }
    int getNormDir() { return norm_dir; }

    private:
    int dim;
    apf::FieldShape* fs;
    apf::EntityShape* es;
    double vol;
    int norm_dir;
  };
  double measureModelEntity(apf::ModelEntity* mdl_ent, apf::Mesh* msh)
  {
    double vol = 0.0;
    int dim = msh->getModelType(mdl_ent);
    apf::MeshEntity* ent;
    auto* it = msh->begin(dim);
    while ((ent = msh->iterate(it))) {
      apf::ModelEntity* to_compare = msh->toModel(ent);
      if (to_compare == mdl_ent) {
        apf::MeshElement* mnt = apf::createMeshElement(msh, ent);
        vol += apf::measure(mnt);
        apf::destroyMeshElement(mnt);
      }
    }
    msh->end(it);
    vol = amsi::comm_sum(vol);
    return vol;
  }

  double measureDisplacedModelEntity(apf::ModelEntity * mdl_ent,
                                     apf::Field * u,
                                     apf::Numbering * numbering)
  {
    double vol = 0.0;
    MeasureDisplaced elvol(u, numbering, 1);
    apf::Mesh* msh = apf::getMesh(u);
    int dim = msh->getModelType(mdl_ent);
    apf::MeshEntity* ent;
    auto* it = msh->begin(dim);
    while ((ent = msh->iterate(it))) {
      apf::ModelEntity* to_compare = msh->toModel(ent);
      if (to_compare == mdl_ent) {
        apf::MeshElement* mnt = apf::createMeshElement(msh, ent);
        elvol.process(mnt);
        vol += elvol.getVol();
        apf::destroyMeshElement(mnt);
      }
    }
    msh->end(it);
    vol = amsi::comm_sum(vol);
    return vol;
  }

  double measureDisplacedMeshEntity(apf::MeshEntity * ent,
                                    apf::Field * u,
                                    apf::Numbering * numbering)
  {
    // todo : derive integration order from field order
    MeasureDisplaced elemental_volume(u, numbering, 1);
    apf::Mesh* msh = apf::getMesh(u);
    apf::MeshElement* mlm = apf::createMeshElement(msh, ent);
    elemental_volume.process(mlm);
    double vol = elemental_volume.getVol();
    apf::destroyMeshElement(mlm);
    return vol;
  }

  double measureDisplacedMeshEntity_greens(apf::MeshEntity * ent,
                                           apf::Field * u,
                                           int norm_dir,
                                           apf::Numbering * numbering)
  {
    // todo : derive integration order from field order
    MeasureDisplacedFromSurf elemental_volume(u, numbering, 1, norm_dir);
    apf::Mesh* msh = apf::getMesh(u);
    apf::MeshElement* mlm = apf::createMeshElement(msh, ent);
    elemental_volume.process(mlm);
    double volume = elemental_volume.getVol();
    apf::destroyMeshElement(mlm);
    return volume;
  }
  double measureDisplacedModelEntity_greens(apf::ModelEntity* mdl_ent,
                                            apf::Field* u)
  {
    double vol = 0.0;
    apf::Mesh* msh = apf::getMesh(u);
    int dim = msh->getModelType(mdl_ent) - 1;
    apf::MeshEntity* ent;
    auto* it = msh->begin(dim);
    while ((ent = msh->iterate(it))) {
      if (msh->isOwned(ent)) {
        apf::ModelEntity* to_compare = msh->toModel(ent);
        if (to_compare == mdl_ent) {
          int sd = side(mdl_ent, msh, ent);
          // assert(sd == 1);
          vol += measureDisplacedMeshEntity_greens(ent, u, sd, nullptr);
        }
      }
    }
    msh->end(it);
    vol = amsi::comm_sum(vol);
    return vol;
  }
}  // namespace amsi