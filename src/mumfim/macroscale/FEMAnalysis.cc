#include "mumfim/macroscale/FEMAnalysis.h"

#include <amsiPETScLAS.h>
#include <apf.h>
#include <petscsnes.h>

#include <iostream>

#include "AnalysisStep.h"
#include "mumfim/exceptions.h"
#include "mumfim/macroscale/ModelTraits.h"
#include "mumfim/macroscale/PetscSNES.h"

// #define MumfimPetscCall(petsc_error_code) if(petsc_error_code) [[unlikely]] {
//  throw mumfim::petsc_error(petsc_error_code); }
namespace mumfim
{
  /**
   * computes the residual vector of the FEM system using the finite element method
   * \param s snes solver object
   * \param solution solution vector to use as an input to calculate the residual
   * \param residual the residual vector computed from the FEM solution
   * \param ctx for this function, the ctx object refers to the FEMAnalysis
   */
  PetscErrorCode FEMAnalysis::CalculateResidual(::SNES s,
                                                Vec solution,
                                                Vec residual,
                                                void * ctx)
  {
    PetscInt iteration;
    SNESGetIterationNumber(s, &iteration);
    // bit hacky...if the last finalized iteration is same as
    // previous
    static int num_calls = 0;
    auto * an = static_cast<FEMAnalysis *>(ctx);
    auto * petsc_las = dynamic_cast<amsi::PetscLAS *>(an->las);
    try
    {
      // in this case, we have called Function another time without
      // checking for convergence. We sent 0 state to tell the
      // microscale that the previous step was not "accepted"
      if (an->iteration > 0)
      {
        an->finalizeIteration(0);
      }
      // Note iteration is not here an actual count of the
      // iterations performed instead it is a test to make sure we
      // call finalize iteration in the correct order between here
      // and the call to check convergence
      ++an->iteration;
      // Given the trial displacement x, compute the residual
      // UpdateDOF (Nonlinear Tissue)
      // 1. Write the new solution into the displacement field
      // las->iter() // we don't need to do this step since we are
      // just using the K/r storage in las not actually computing
      // residuals
      const double * sol;
      VecGetArrayRead(solution, &sol);
      an->analysis_step_->UpdateDOFs(sol);
      VecRestoreArrayRead(solution, &sol);
      an->analysis_step_->ApplyBC_Dirichlet();
      an->analysis_step_->RenumberDOFs();
      int gbl, lcl, off;
      an->analysis_step_->GetDOFInfo(gbl, lcl, off);
      an->las->Reinitialize(gbl, lcl, off);
      an->las->Zero();
      // assembles into Mat/vec
      an->analysis_step_->Assemble(an->las);
      VecAssemblyBegin(petsc_las->GetVector());
      VecAssemblyEnd(petsc_las->GetVector());
      MatAssemblyBegin(petsc_las->GetMatrix(), MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(petsc_las->GetMatrix(), MAT_FINAL_ASSEMBLY);
      an->analysis_step_->iter();
      VecCopy(petsc_las->GetVector(), residual);
      // an->analysis_step_->AcceptDOFs();
    }
    catch (mumfim_error & e)
    {
      std::cerr << "Took bad step. Writing debug output\n";
      std::cerr << e.what() << "\n";
      an->stp = -2;
      an->checkpoint();
      PetscViewer viewer;
      // write residual/stiffness to file
      auto stiffness_matrix_file =
          std::string(amsi::fs->getResultsDir() + "/petsc_vector_states.mtx");
      MumfimPetscCall(PetscViewerBinaryOpen(AMSI_COMM_SCALE,
                                            stiffness_matrix_file.c_str(),
                                            FILE_MODE_WRITE, &viewer));
      // tangent stiffness matrix
      MumfimPetscCall(MatView(petsc_las->GetMatrix(), viewer));
      // residual vector
      MumfimPetscCall(VecView(petsc_las->GetVector(), viewer));
      MumfimPetscCall(SNESSetFunctionDomainError(s));
    }
    return 0;
  }

  /**
   * computes the Jacobian Matrix of the FEM system using the finite element
   * method.
   * \warning this function does not actually compute the jacobian, but copies
   * the jacobian that is constructed when the residual is computed
   * \param s snes solver object
   * \param solution solution vector to use as an input to calculate the
   * residual
   * \param residual the residual vector computed from the FEM solution
   * \param ctx for this function, the ctx object refers to the FEMAnalysis
   */
  PetscErrorCode FEMAnalysis::CalculateJacobian(::SNES /* unused (snes) */,
                                                Vec /* unused (solution) */,
                                                Mat Amat,
                                                Mat /* unused (PMat) */,
                                                void * ctx)
  {
    auto * an = static_cast<FEMAnalysis *>(ctx);
    auto * petsc_las = dynamic_cast<amsi::PetscLAS *>(an->las);
    MumfimPetscCall(MatCopy(petsc_las->GetMatrix(), Amat, SAME_NONZERO_PATTERN));
    MumfimPetscCall(MatScale(Amat, -1));
    return 0;
  }

  PetscErrorCode FEMAnalysis::CheckConverged(::SNES snes,
                                             PetscInt it,
                                             PetscReal xnorm,
                                             PetscReal gnorm,
                                             PetscReal f,
                                             SNESConvergedReason * reason,
                                             void * ctx)
  {
    auto * an = static_cast<FEMAnalysis *>(ctx);
    auto error = SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, ctx);
    bool converged = (reason != nullptr && *reason != SNES_CONVERGED_ITERATING);
    int accepted = converged ? 1 : -1;
    // For MultiscaleAnalysis this informs microscale if the step is
    // done the microscale knows that a value of 0 means step is
    // accepted but not converged
    an->finalizeIteration(accepted);
    an->analysis_step_->AcceptDOFs();
    // HACK set this to zero so we don't finalize iteration
    // in form function
    an->iteration = 0;
    return error;
  }

  FEMAnalysis::FEMAnalysis(apf::Mesh * mesh,
                                 std::unique_ptr<const mt::CategoryNode> cs,
                                 MPI_Comm c,
                                 const amsi::Analysis & amsi_analysis)
      : cm(c)
      , analysis_case(std::move(cs))
      , mesh(mesh)
      , t(0.0)
      , dt(0.0)
      , stp(0)
      , mx_stp(1)
      , analysis_step_(nullptr)
      , las(new amsi::PetscLAS(0, 0))
      , completed(false)
      , state_fn()
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
    mx_stp = (*timesteps_trait)();
    dt = (double)1.0 / (double)mx_stp;
    // output params
#ifdef LOGRUN
    int rnk = -1;
    std::stringstream cnvrt;
    MPI_Comm_rank(cm, &rnk);
    cnvrt << rnk;
    state_fn =
        amsi::fs->getResultsDir() + "/tissue_state." + cnvrt.str() + ".log";
    state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << "STEP, ITER,   T, DESC\n"
                     << "   0,    0, 0.0, init\n";
#endif
  }
  FEMAnalysis::~FEMAnalysis()
  {
    delete analysis_step_;
    delete las;
#ifdef LOGRUN
    amsi::deleteLog(state);
#endif
  }
  void FEMAnalysis::run()
  {
    try
    {
      analysis_step_->preRun();  // calls updateMicro for multiscale analysis
      analysis_step_->recoverSecondaryVariables(stp);
      checkpoint();
      // write the initial state of everything
      t += dt;
      analysis_step_->setSimulationTime(t);
      // analysis_step_->getUField());
      analysis_step_->computeInitGuess(las);
      completed = false;
      SNES snes(cm);
      MumfimPetscCall(SNESSetFromOptions(snes));
      auto * petsc_las = dynamic_cast<amsi::PetscLAS *>(las);
      if (petsc_las == nullptr)
      {
        std::cerr << "Current solver only works with petsc backend!\n";
        std::exit(1);
      }
      // set up the SNES functions
      Mat AMat;
      MumfimPetscCall(
          MatDuplicate(petsc_las->GetMatrix(), MAT_DO_NOT_COPY_VALUES, &AMat));
      MumfimPetscCall(SNESSetFunction(snes, nullptr, CalculateResidual,
                                      static_cast<void *>(this)));
      MumfimPetscCall(SNESSetJacobian(snes, AMat, AMat, CalculateJacobian,
                                      static_cast<void *>(this)));
      MumfimPetscCall(SNESSetConvergenceTest(
          snes, CheckConverged, static_cast<void *>(this), nullptr));

      while (!completed)
      {
#ifdef LOGRUN
        amsi::log(state) << stp << ", " << MPI_Wtime() << ", "
                         << "start_step" << std::endl;
#endif
        if (!PCU_Comm_Self()) std::cout << "Load step = " << stp << std::endl;
        // TODO wrap SNES in RAII class so create/destroy is exception safe
        MumfimPetscCall(
            SNESSolve(snes, nullptr, petsc_las->GetSolutionVector()));
        const auto converged = std::invoke(
            [&snes]()
            {
              SNESConvergedReason converged;
              MumfimPetscCall(SNESGetConvergedReason(snes, &converged));
              return converged;
            });
        const auto iterations = std::invoke(
            [&snes]()
            {
              PetscInt iteration;
              MumfimPetscCall(SNESGetIterationNumber(snes, &iteration));
              return iteration;
            });
        // analysis diverged
        if (converged < 0)
        {
          completed = true;
          // finalizeStep(); // shouldn't call this here since it's called at
          // the end
          std::stringstream ss;
          ss << "Step " << stp << "failed to converge\n";
          ss << SNESConvergedReasons[converged];
          ss << "Number of nonlinear iterations = " << iterations << "\n";
          throw mumfim_error(ss.str());
        }
        if (stp >= mx_stp - 1)
        {
          completed = true;
          std::cout << "Final load step converged. Case complete." << std::endl;
        }
        std::cout << "checkpointing (macro)" << std::endl;
        std::cout << "Rewriting at end of load step to include orientation data"
                  << std::endl;
        analysis_step_->recoverSecondaryVariables(stp);
        checkpoint();
        stp++;
        t += dt;
        analysis_step_->setSimulationTime(t);
        // Warning! this function has a potentially blocking MPI CALL!
        finalizeStep();
      }
    }
    catch (mumfim_error & e)
    {
      stp = -4;
      // Write out the failed state
      // analysis_step_->recoverSecondaryVariables(stp);
      checkpoint();
      PetscViewer viewer;
      auto * petsc_las = dynamic_cast<amsi::PetscLAS *>(this->las);
      if (petsc_las)
      {
        // write residual/stiffness to file
        auto stiffness_matrix_file =
            std::string(amsi::fs->getResultsDir() + "/petsc_vector_states.mtx");
        MumfimPetscCall(PetscViewerBinaryOpen(AMSI_COMM_SCALE,
                                              stiffness_matrix_file.c_str(),
                                              FILE_MODE_WRITE, &viewer));
        // tangent stiffness matrix
        MumfimPetscCall(MatView(petsc_las->GetMatrix(), viewer));
        // residual vector
        MumfimPetscCall(VecView(petsc_las->GetVector(), viewer));
      }
      throw;
    }
  }
  void FEMAnalysis::finalizeStep(){};
  void FEMAnalysis::finalizeIteration(int){};
  void FEMAnalysis::checkpoint()
  {
#ifdef LOGRUN
    std::ofstream st_fs(state_fn.c_str(), std::ios::out | std::ios::app);
    amsi::flushToStream(state, st_fs);
#endif
    std::stringstream cnvrt;
    cnvrt << "msh_stp_" << stp;
    apf::writeVtkFiles(
        std::string(amsi::fs->getResultsDir() + "/" + cnvrt.str()).c_str(),
        analysis_step_->getMesh());
  }
}  // namespace mumfim
