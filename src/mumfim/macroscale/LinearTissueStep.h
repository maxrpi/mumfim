#ifndef MUMFIM_LINEARTISSUE_H_
#define MUMFIM_LINEARTISSUE_H_
#include <memory>

#include "amsiFEA.h"

namespace mumfim
{
  class LinearTissueStep : public amsi::FEAStep
  {
    protected:
    std::map<int, std::unique_ptr<amsi::ElementalSystem> > constitutives;
    public:
    LinearTissueStep(apf::Mesh * mesh,
                    const mt::CategoryNode & _analysis_case,
                    MPI_Comm cm);
    virtual void UpdateDOFs(const double * sol) override;
    virtual void Assemble(amsi::LAS * las) override;
    [[nodiscard]] apf::Field * getUField() final { return apf_primary_field; }

    [[nodiscard]] amsi::ElementalSystem * getIntegrator(
        apf::MeshEntity * mesh_entity,
        int integration_point) override;
  };
}  // namespace mumfim
#endif
