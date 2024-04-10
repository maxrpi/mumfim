#ifndef AMSI_ANALYSIS_H_
#define AMSI_ANALYSIS_H_
#include <amsiAnalysisConfig.h>
#include <amsiPETScLAS.h>
#include <amsiNonlinearAnalysis.h>
#include <amsiConfigurationOptions.h>
#include <amsiMPI.h>
#include <amsiFileSys.h>
#include <amsiUtil.h>
namespace amsi {
  class Analysis {
    class PetscScopeGuard {
      public:
      PetscScopeGuard() = default;
      [[nodiscard]] PetscScopeGuard(int argc, char** argv,
                                    const std::string& petsc_options_file,
                                    MPI_Comm cm);
      ~PetscScopeGuard();
      PetscScopeGuard(const PetscScopeGuard&) = delete;
      PetscScopeGuard(PetscScopeGuard&& old) noexcept
      {
        std::swap(initialized_here_, old.initialized_here_);
      }
      PetscScopeGuard& operator=(const PetscScopeGuard&) = delete;
      PetscScopeGuard& operator=(PetscScopeGuard&& old) noexcept
      {
        std::swap(initialized_here_, old.initialized_here_);
        return *this;
      }

      private:
      bool initialized_here_{false};
    };

    public:
    Analysis(const AnalysisOptions& options, int argc, char** argv, MPI_Comm cm,
             const MPI& mpi)
        : results_directory_(
              std::make_unique<FileSystemInfo>(options.results_directory))
        , cm_(cm), mpi_(mpi)
    {
      // workaround until we get rid of the global variables!
      SetGlobalFilesystem(results_directory_.get());
      if (options.use_petsc) {
        petsc_scope_guard_ =
            PetscScopeGuard{argc, argv, options.petsc_options_file, cm};
      }
    }
    [[nodiscard]] const std::string& getResultsDir() const noexcept
    {
      return results_directory_->getResultsDir();
    }
    [[nodiscard]] const MPIComm& getComm() const noexcept { return cm_; }

    private:
    [[maybe_unused]] PetscScopeGuard petsc_scope_guard_;
    std::unique_ptr<FileSystemInfo> results_directory_;
    MPIComm cm_;
    const MPI& mpi_;
  };
}  // namespace amsi
#endif
