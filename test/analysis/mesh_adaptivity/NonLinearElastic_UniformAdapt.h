#ifndef NONLINEARELASTIC_UNIFORMADAPT_H_
#define NONLINEARELASTIC_UNIFORMADAPT_H_
#include <apfSIM.h>
#include <NonLinElasticity.h>
#include <cassert>
namespace amsi
{
  class UniformAdapt : public NonLinElasticity
  {
  public:
    UniformAdapt(pGModel imdl,
                 pParMesh imsh,
                 pACase pd,
                 MPI_Comm cm = AMSI_COMM_SCALE)
      : FEA(cm)
      , NonLinElasticity(imdl,imsh,pd,cm)
    {
      mesh_size_field = apf::createSIMLagrangeField(apf_mesh,"mesh_size",apf::SCALAR,1);
    }
    virtual void Adapt()
    {
      assert(mesh_size_field);
      should_adapt = true;
      if(should_adapt)
      {
        apf::MeshEntity * vert = NULL;
        for(apf::MeshIterator * it = apf_mesh->begin(0); (vert = apf_mesh->iterate(it));)
          apf::setScalar(mesh_size_field,vert,0,0.1);
        apfSimFEA::Adapt();
        numbered = false;
        should_adapt = false;
      }
    }
  };
}
#endif
