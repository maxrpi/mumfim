#ifndef APF_FIELD_OP_H_
#define APF_FIELD_OP_H_
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <vector>
namespace amsi
{
  class FieldOp
  {
  public:
    virtual bool inEntity(apf::MeshEntity *) {return true;}
    virtual void outEntity() {}
    virtual void atNode(int) {}
    void apply(apf::Field * f)
    {
      apf::Mesh * msh = apf::getMesh(f);
      apf::FieldShape * s = apf::getShape(f);
      for(int ii = 0; ii < 4; ++ii)
      {
        if(!s->hasNodesIn(ii))
          continue;
        apf::MeshIterator * it = msh->begin(ii);
        apf::MeshEntity * e = NULL;
        while((e = msh->iterate(it)))
        {
          int n = s->countNodesOn(msh->getType(e));
          if(n==0)
            continue;
          if(!this->inEntity(e))
            continue;
          for(int jj = 0; jj < n; ++jj)
            this->atNode(jj);
          this->outEntity();
        }
        msh->end(it);
      }
    }
    // field op apply over structures of apf::MeshEntity* e.g. vector<apf::MeshEntity*>
    void apply(std::vector<apf::MeshEntity *>::iterator ent_bgn,
               std::vector<apf::MeshEntity *>::iterator ent_end,
               apf::Field * f)
    {
      apf::Mesh * msh = apf::getMesh(f);
      apf::FieldShape * s = apf::getShape(f);
      for (int ii = 0; ii < 4; ++ii)
      {
        if (!s->hasNodesIn(ii)) continue;
        for (auto e = ent_bgn; e != ent_end; ++e)
        {
          int n = s->countNodesOn(msh->getType(*e));
          if (n == 0) continue;
          if (!this->inEntity(*e)) continue;
          for (int jj = 0; jj < n; ++jj)
            this->atNode(jj);
          this->outEntity();
        }
      }
    }
  };
}
#endif
