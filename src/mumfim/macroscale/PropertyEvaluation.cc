#include "mumfim/macroscale/PropertyEvaluation.h"
#include "StepperFactory.h"

#include <amsiPETScLAS.h>
#include <apf.h>
#include <petscsnes.h>
#include <petscvec.h>

#include <iostream>

#include "amsiFEA.h"
#include "mumfim/exceptions.h"
#include "mumfim/macroscale/ModelTraits.h"
#include "mumfim/macroscale/PetscSNES.h"

// #define MumfimPetscCall(petsc_error_code) if(petsc_error_code) [[unlikely]] {
//  throw mumfim::petsc_error(petsc_error_code); }
namespace mumfim
{

  PropertyEvaluation::PropertyEvaluation(apf::Mesh * mesh,
                                 std::unique_ptr<const mt::CategoryNode> cs,
                                 std::string ktf,
                                 MPI_Comm c,
                                 const amsi::Analysis & amsi_analysis)
      : cm(c)
      , analysis_case(std::move(cs))
      , mesh(mesh)
      , analysis_step_(nullptr)
      , kappa_tag_filename(ktf)
  {
    // util data
    const auto * problem_definition =
        mt::GetPrimaryCategoryByType(analysis_case.get(), "problem definition");
    const auto * solution_strategy =
        mt::GetPrimaryCategoryByType(analysis_case.get(), "solution strategy");
    if (problem_definition == nullptr || solution_strategy == nullptr ||
        problem_definition->GetType() != "macro" ||
        solution_strategy->GetType() != "macro")
    {
      std::cerr << "Analysis case should have  \"problem definition\" and "
                   "\"solution strategy\" of the \"macro\" analysis type.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    // analysis params
    const auto * timesteps_trait = mt::GetCategoryModelTraitByType<mt::IntMT>(
        solution_strategy, "num timesteps");
    if (timesteps_trait == nullptr)
    {
      std::cerr << R"("solution strategy" must have "num timesteps" trait)";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    analysis_step_ = createStepper(mesh, *analysis_case, cm, kappa_tag_filename);
  }
  PropertyEvaluation::~PropertyEvaluation()
  {
    delete analysis_step_;
  }
  void PropertyEvaluation::run()
  {
    try
    {
      analysis_step_->Assemble(nullptr);
      // analysis_step_->Solve();
    }
    catch (mumfim_error & e)
    {
      std::cout << "Error in run()\n";
    }
    
  }
}  // namespace mumfim
