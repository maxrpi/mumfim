#include "amsiAnalysis.h"
namespace amsi {
  Analysis::PetscScopeGuard::PetscScopeGuard(
      int argc, char **argv, const std::string &petsc_options_file, MPI_Comm cm)
  {
#ifdef PETSC
    PetscBool initialized = PETSC_FALSE;
    // PetscInitialized(&initialized);
    //PETSC_COMM_WORLD = cm;
    if (initialized == PETSC_FALSE) {
      initialized_here_ = true;
      PetscInitialize(&argc, &argv, petsc_options_file.c_str(), PETSC_NULL);
    }
#endif
  }
  Analysis::PetscScopeGuard::~PetscScopeGuard()
  {
#ifdef PETSC
    if (initialized_here_) {
      PetscFinalize();
    }
#endif
  }
}  // namespace amsi
