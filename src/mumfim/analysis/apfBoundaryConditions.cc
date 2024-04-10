#include "apfBoundaryConditions.h"
#include "apfShape.h"
#include "apfNumbering.h"
#include <iostream>
namespace amsi
{
  double getDirichletValue(BCQuery * qry,
                           apf::Mesh * msh,
                           apf::MeshEntity * ent,
                           int nd,
                           int cmp,
                           double t)
  {
    double val = 0.0;
    bool cnst = qry->isConst(cmp);
    bool sptl = qry->isSpaceExpr(cmp);
    bool time = qry->isTimeExpr(cmp);
    if(cnst)
      val = qry->getValue(cmp);
    else
    {
      if(sptl)
      {
        // get parametric coordinate of nd
        apf::Vector3 pt;
        msh->getPoint(ent,nd,pt);
        // map local coord to global coord
        apf::Vector3 xyz;
        apf::MeshElement * melmt =apf::createMeshElement(msh,ent);
        apf::mapLocalToGlobal(melmt,pt,xyz);
        apf::destroyMeshElement(melmt);
        // check if also time expr
        if(time)
          val = qry->getValue(cmp,t,xyz[0],xyz[1],xyz[2]);
        else
          val = qry->getValue(cmp,xyz[0],xyz[1],xyz[2]);
      }
      else if(time)
        val = qry->getValue(cmp,t);
      else
        std::cerr << "Error: Non constant dirichlet boundary condition with "
                  << " no temporal or spatial dependencies." << std::endl;
    }
    return val;
  }
}
