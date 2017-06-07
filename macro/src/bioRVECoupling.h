#ifndef BIO_RVE_COUPLING_H_
#define BIO_RVE_COUPLING_H_
#include <bioMicroFOMultiscale.h>
namespace bio
{
  enum MicroscaleType
  {
    NONE = 0,
    FIBER_ONLY = 1,
    FIBER_MATRIX = 2,
    MICROSCALE_TYPE_COUNT = 3
  };
  class RVECoupling
  {
  private:
    apf::Mesh * msh;
    apf::Field * rst_fst;
    apf::Field * crt_rve;
    apf::Field * prv_rve;
    int fld_ord;
    int rve_cnt;
    int rve_rgns;
    std::vector<micro_fo_result> rsts;
  public:
  RVECoupling(apf::Mesh * m, int o)
      : msh(m)
      , rst_fst(NULL)
      , crt_rve(NULL)
      , prv_rve(NULL)
      , fld_ord(o)
      , rve_cnt(0)
      , rve_rgns(0)
      , rsts()
    {
      rst_fst = apf::createIPField(msh,"micro_fo_rve_offset",apf::SCALAR,1);
      crt_rve = apf::createIPField(msh,"micro_fo_current_rve",apf::SCALAR,1);
      prv_rve = apf::createIPField(msh,"micro_fo_previous_rve",apf::SCALAR,1);
    }
    ~RVECoupling()
    {
      apf::destroyField(prv_rve);
      apf::destroyField(crt_rve);
      apf::destroyField(rst_fst);
    }
    micro_fo_result * getRVEResult(apf::MeshEntity * me, int ip)
    {
      int fst = apf::getScalar(rst_fst,me,ip);
      return &rsts[fst];
    }
    int countRVEsOn(apf::MeshEntity * me)
    {
      int cnt = 0;
      if(apf::getScalar(crt_rve,me,0) == FIBER_ONLY)
      {
        apf::MeshElement * mlm = apf::createMeshElement(msh,me);
        cnt = apf::countIntPoints(mlm,fld_ord);
        apf::destroyMeshElement(mlm);
      }
      return cnt;
    }
  };
}
#endif
