#ifndef AMSI_AMSIBOUNDARYCONDITIONS_H
#define AMSI_AMSIBOUNDARYCONDITIONS_H
#include <apfNumbering.h>
#include <model_traits/AssociatedModelTraits.h>
#include "amsiLAS.h"
#include "amsiNeumannIntegratorsMT.h"
namespace amsi {
  using bc_path_type =
      std::vector<std::pair<std::string, std::vector<std::string>>>;
  struct DirichletBCEntry {
    std::vector<std::string> categories;
    std::string mt_name;
  };
  struct NeumannBCEntry {
    std::vector<std::string> categories;
    std::string mt_name;
    NeumannBCType mt_type;
  };
  void applyNeumannBCs(LAS* las, apf::Numbering* nm,
                       const mt::AssociatedModelTraits<mt::DimIdGeometry>& bcs,
                       const std::vector<NeumannBCEntry>& bc_paths, double t);
  int applyDirichletBCs(apf::Numbering* nm, apf::Field* delta_field,
                        const mt::AssociatedModelTraits<mt::DimIdGeometry>& bcs,
                        const std::vector<DirichletBCEntry>& bc_paths,
                        double t);
}  // namespace amsi
#endif  // AMSI_AMSIBOUNDARYCONDITIONS_H
