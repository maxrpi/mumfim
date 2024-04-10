#ifndef SPRADAPT_H_
#define SPRADAPT_H_

#include <apfSIM.h>
#include <NonLinElasticity.h>

#include <spr.h>

namespace amsi {

  class SPRAdapt : public Analysis::NonLinElasticity
  {
  private:
    apf::Field * stress_ip_field;
    
  public:
    SPRAdapt(MPI_Comm comm,
	     pGModel in_model,
	     pParMesh in_mesh) :
    FEA(comm,"uniform_adapt"),
    NonLinElasticity(comm,in_model,in_mesh)
    {
      mesh_size_field = apf::createSIMLagrangeField(apf_mesh,"mesh_size",apf::SCALAR,1);
      stress_ip_field = apf::createIPField(apf_mesh,"stress_ip_field",apf::MATRIX,1);
    }
    
    virtual void Adapt()
    {
      assert(mesh_size_field);

      should_adapt = true;
      if(should_adapt)
      {
	mesh_size_field = spr::getSPRSizeField(stress_ip_field,0.1);
	apfSimFEA::Adapt();
	numbered = false;
	should_adapt = false;
	WriteMesh("post_spr_adapt_mapped");
      }
    }
  };
  
}

#endif
