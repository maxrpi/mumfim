#ifndef APF_MESH_ITERATOR_H_
#define APF_MESH_ITERATOR_H_
#include <apf.h>
#include <apfMesh.h>
namespace amsi
{
  class apfMeshIterator
  {
  protected:
    apf::Mesh * msh; // unowned
    apf::MeshIterator * itr;
    apf::MeshEntity * crt;
    int ent_dim;
    apfMeshIterator()
      : msh(NULL)
      , itr(NULL)
      , crt(NULL)
      , ent_dim(-1)
    { }
  public:
    apfMeshIterator(apf::Mesh * m, int ed)
      : msh(m)
      , itr(NULL)
      , crt(NULL)
      , ent_dim(ed)
    {
      itr = msh->begin(ed);
      operator++();
    }
    apfMeshIterator(const apfMeshIterator & o)
      : apfMeshIterator(o.msh,o.ent_dim) //c++11
    { }
    virtual ~apfMeshIterator()
    {
      msh->end(itr);
    }
    virtual void operator++()
    {
      crt = msh->iterate(itr);
    }
    virtual bool operator==(const apfMeshIterator & o) const
    {
      return crt == o.crt && msh == o.msh;
    }
    virtual bool operator!=(const apfMeshIterator & o) const { return !operator==(o); }
    virtual apf::MeshEntity * operator*() const { return crt; }
  };
  class apfEndIterator : public apfMeshIterator
  {
  public:
    apfEndIterator(apf::Mesh * m)
      : apfMeshIterator()
    {
      msh = m;
    }
    virtual ~apfEndIterator() { }
    virtual void operator++() { }
  };
}
#endif
