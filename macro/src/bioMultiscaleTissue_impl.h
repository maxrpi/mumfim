#include "amsiDeformation.h"
namespace bio
{
  /*
    template <typename O>
    void MultiscaleTissue::updateRVEDeletion(O o, bool all)
    {
    int iid = 0;
    for(auto me = rve_ents.begin(); me != rve_ents.end(); me++)
    {
    apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,*me);
    int num_gauss_pts = apf::countIntPoints(mlm,getOrder(mlm));
    for(int ii = 0 ; ii < num_gauss_pts; ii++)
    {
    int crt = apf::getScalar(crt_rve,*me,ii);
    int prv = apf::getScalar(prv_rve,*me,ii);
    if((crt == NONE && prv != NONE) || all)
    {
    o = iid;
    ++o;
    // assumed that removing one RVE from an entity removed that entity from the list of entities containing RVEs, not necessarily the case when using multiply IP per ent
    //me = rve_ents.erase(me);
    //me--;
    }
    iid++;
    }
    apf::destroyMeshElement(mlm);
    }
    }
  */
  template <typename O1, typename O2, typename O3>
    void MultiscaleTissue::serializeNewRVEData(O1 new_hdrs, O2 new_prms, O3 new_data, bool all)
  {
    apf::MeshEntity * rgn = NULL;
    apf::MeshIterator * it = NULL;
    for(it = apf_mesh->begin(3); (rgn = apf_mesh->iterate(it));)
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,rgn);
      int ng = apf::countIntPoints(mlm,getOrder(mlm));
      for(int ip = 0; ip < ng; ++ip)
      {
        int crt = apf::getScalar(crt_rve,rgn,ip);
        int prv = apf::getScalar(prv_rve,rgn,ip);
        if((crt == FIBER_ONLY && prv != FIBER_ONLY) || all) // if the RVE is new
        {
          micro_fo_header hdr;
          micro_fo_params prm;
          micro_fo_init_data dat;
          getInternalRVEData(rgn,hdr,prm,dat);
          getExternalRVEData(rgn,hdr,prm);
          hdr.data[GAUSS_ID] = ip;
          *new_hdrs++ = hdr;
          *new_prms++ = prm;
          *new_data++ = dat;
        }
      }
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
  }
  // TODO switchover to using deformation gradient means that we
  //  oneed to know which integration point in the element we are
  //  talking about when we serialize the RVE data, previously
  //  micro already knew this from initialization so it wasn't
  //  important at macro
  template <typename O>
    void MultiscaleTissue::serializeRVEData(O o)
  {
    apf::MeshEntity * rgn = NULL;
    for(auto * it = apf_mesh->begin(3); (rgn = apf_mesh->iterate(it)); )
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh,rgn);
      apf::Element * e = apf::createElement(apf_primary_field,mlm);
      int ng = apf::countIntPoints(mlm,getOrder(mlm));
      for(int ip = 0; ip < ng; ++ip)
      {
        int crt = apf::getScalar(crt_rve,rgn,ip);
        if(crt == FIBER_ONLY)
        {
          apf::Matrix3x3 F;
          apf::Vector3 p;
          apf::getIntPoint(mlm,1,ip,p);
          amsi::deformationGradient(e,p,F);
          micro_fo_data data;
          for(int ii = 0; ii < 3; ++ii)
            for(int jj = 0; jj < 3; ++jj)
              data.data[ii*3 + jj] = F[ii][jj];
          *o++ = data;
        }
      }
      apf::destroyElement(e);
      apf::destroyMeshElement(mlm);
    }
  }
  template <typename I>
    int getRVEDirectoryIndex(I tp_bgn, I tp_end, apf::ModelEntity * ent)
  {
    int ii = -1;
    pGEntity rgn = reinterpret_cast<pGEntity>(ent);
    pAttribute mdl = GEN_attrib(rgn,"material model");
    pAttribute sm = Attribute_childByType(mdl,"multiscale model");
    if(sm)
    {
      pAttributeString dir  = (pAttributeString)Attribute_childByType(sm,"directory");
      pAttributeString prfx = (pAttributeString)Attribute_childByType(sm,"prefix");
      char * dir_str = AttributeString_value(dir);
      char * tp_str  = AttributeString_value(prfx);
      std::string tp(std::string(dir_str) + std::string("/") + std::string(tp_str));
      auto fnd = std::find(tp_bgn,tp_end,tp);
      if(fnd != tp_end)
        ii = std::distance(tp_bgn,fnd);
      Sim_deleteString(dir_str);
      Sim_deleteString(tp_str);
    }
    return ii;
  }
}
