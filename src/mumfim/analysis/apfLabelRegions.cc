#include <apf.h>
#include <apfMesh.h>
namespace amsi
{
  void applyUniqueRegionTags(apf::Mesh* mesh)
  {
    auto * fld = apf::createIPField(mesh,"region_id",apf::SCALAR,1);
    apf::MeshEntity* ent;
    auto* it = mesh->begin(3);
    while((ent = mesh->iterate(it)))
    {
      auto* model_region = mesh->toModel(ent);
      apf::setScalar(fld,ent,0, mesh->getModelTag(model_region));

    }
  }
}
