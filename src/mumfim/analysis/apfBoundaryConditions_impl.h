#ifndef APF_BOUNDARY_CONDITIONS_IMPL_H_
#define APF_BOUNDARY_CONDITIONS_IMPL_H_
#include <cassert>
#include "amsiBoundaryConditionQuery.h"
#include "amsiNeumannIntegrators.h"
#include "apfNumbering.h"
#include "apfShape.h"
namespace amsi
{
  template <typename I>
    int applyDirichletBC(apf::Numbering * nm,
                         const I& begin,
                         const I& end,
                         BCQuery * qry,
                         double t,
                         bc_set_map & already_set_map,
                         apf::Field* deltaField)
  {
    int fxd = 0;
    apf::Field * fld = apf::getField(nm);
    apf::Mesh * msh = apf::getMesh(fld);
    apf::FieldShape * fs = apf::getShape(fld);
    int cmps = apf::countComponents(fld);
    double * vls = new double[cmps];
    double * delta_vls = new double[cmps];
    assert(qry->numComps() == cmps);
    for(I it = begin; it != end; ++it)
    {
      // TODO : would prefer to implement a wrapper on the sim iterator...
      apf::MeshEntity * ent = reinterpret_cast<apf::MeshEntity*>(*it);
      int nds = fs->countNodesOn(msh->getType(ent));
      for(int ii = 0; ii < nds; ii++)
      {
        apf::getComponents(fld,ent,ii,vls);
        if(deltaField)
          apf::getComponents(deltaField,ent,ii,delta_vls);
        for(int jj = 0; jj < cmps; jj++)
        {
          if(qry->isFixed(jj))
          {
            auto has_ent = already_set_map.find(ent);
            // if the entity doesn't have any bc's set yet,
            // or the specific component of the given node
            // isn't set yet, then apply the boundary condition
            if(has_ent == already_set_map.end() || 
               !has_ent->second[ii*cmps+jj])
            {
              auto current_loc = already_set_map.emplace(std::make_pair(ent,
                                                 std::vector<bool>(cmps*nds, false)));
              // set the current dof of the current entity to a alredy written state
              current_loc.first->second[ii*cmps+jj] = true;
              apf::fix(nm,ent,ii,jj,true);
              double dirichlet_vl = getDirichletValue(qry,msh,ent,ii,jj,t); 
              // the vls field must already be in the updated state here...
              // how is this possible?
              delta_vls[jj] = dirichlet_vl-vls[jj];
              //delta_vls[jj] = vls[jj];
              vls[jj] = dirichlet_vl;
              fxd++;
            }
          }
        }
        // if the user supplied a delta field such as the delta_U field
        // in a solid mechanics problem set this value as well
        // we should only really set the components that lie of the fixe surfaces
        if(deltaField)
          apf::setComponents(deltaField,ent,ii,delta_vls);
        apf::setComponents(fld,ent,ii,vls);
      }
    }
    delete [] vls;
    delete [] delta_vls;
    return fxd;
  }
  template <typename I>
    void applyNeumannBC(LAS * las,
                        apf::Numbering * nm,
                        const I& bgn,
                        const I& nd,
                        NeumannIntegrator * i,
                        double)
  {
    apf::Field * fld = apf::getField(nm);
    for(I it = bgn; it != nd; ++it)
    {
      apf::MeshEntity * ent = reinterpret_cast<apf::MeshEntity*>(*it);
      apf::MeshElement * mnt = apf::createMeshElement(apf::getMesh(fld),ent);
      i->process(mnt);
      apf::NewArray<int> dofs;
      apf::getElementNumbers(nm,ent,dofs);
      las->AddToVector(i->getnedofs(),&dofs[0],&i->getFe()[0]);
      apf::destroyMeshElement(mnt);
    }
  }
}
#endif
