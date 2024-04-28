#include "StepperFactory.h"
#include "LinearHeatConductionStep.h"
#include "LinearTissueStep.h"
#include "NonlinearTissueStep.h"
#include "MultiscaleTissueStep.h"

namespace mumfim
{
 AnalysisStep * createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      MPI_Comm com
  )
  {
    AnalysisStep *stepper = nullptr;

    const auto * problem_definition =
        mt::GetPrimaryCategoryByType(&analysis_case, "problem definition");
    if (problem_definition == nullptr || problem_definition->GetType() != "macro")
    {
      std::cerr << "Analysis case should have  \"problem definition\" and "
                   "\"macro\" analysis type.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    int problem_type_index = 22; // set default to nonlinear tissue for now,
                                // to not break old problems
    const auto * problem_type_trait = mt::GetCategoryModelTraitByType<mt::IntMT>
          (problem_definition, "problem type");
    if (problem_type_trait == nullptr)
    {
      std::cerr << R"("problem definition" must have "problem type" trait. Defaulting to NonlinearTissue)"<<std::endl;
    }
    else{
      problem_type_index = (*problem_type_trait)();
    }

    switch(problem_type_index){
      case(11):
        stepper = new LinearHeatConductionStep(mesh, analysis_case, com);
        break;
      case(21):
        stepper = new LinearTissueStep(mesh, analysis_case, com);
        break;
      case(22):
        stepper = new NonlinearTissueStep(mesh, analysis_case, com);
        break;
      default:
        std::cerr << "\"problem type\" trait of " << problem_type_index <<
            "is not implemented. Exiting.";
        MPI_Abort(AMSI_COMM_WORLD, 1);
    };
    return stepper;
  }

  AnalysisStep * createStepper(
      apf::Mesh * mesh,
      const mt::CategoryNode & analysis_case,
      const amsi::Multiscale & amsi_multiscale,
      MPI_Comm com
  )
  {
    AnalysisStep *stepper = nullptr;

    const auto * problem_definition =
        mt::GetPrimaryCategoryByType(&analysis_case, "problem definition");
    if (problem_definition == nullptr)
    {
      std::cerr << "Analysis case should have  \"problem definition\".\n ";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * problem_type = mt::GetCategoryModelTraitByType<mt::IntMT>(
        problem_definition, "problem type");
    if (problem_type == nullptr)
    {
      std::cerr << R"("problem definition" must have "problem type" trait)";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    int problem_type_index = (int) (*problem_type)();

    switch(problem_type_index){
      case(3):
        stepper = new MultiscaleTissueStep(mesh, analysis_case, com, amsi_multiscale);
        break;
      default:
        std::cerr << "\"problem type\" trait of " << problem_type_index <<
            "is not implemented for multiscale. Exiting.";
        MPI_Abort(AMSI_COMM_WORLD, 1);
    };
    return stepper;
  }
}