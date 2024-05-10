#include "Solvers.h"
#include "amsiNonlinearAnalysis.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <amsiVerbosity.h>
namespace amsi
{

  void LinearSolver(FEAStep * fem,LAS * las)
  {
    int global_dof_count, local_dof_count, first_local_dof;
    fem->ApplyBC_Dirichlet();
    ///// configure linear system
    fem->RenumberDOFs();
    fem->GetDOFInfo(global_dof_count,local_dof_count,first_local_dof);
    las->Reinitialize(local_dof_count,global_dof_count,first_local_dof);
    las->Zero();
    ///// assemble linear system
    fem->Assemble(las);
    ///// solve linear system
    las->solve();

    /*
    std::cout << "Vector" << std::endl;
    las->PrintVector(std::cout);
    std::cout << "Matrix" << std::endl;
    las->PrintMatrix(std::cout);
    std::cout << "Solution" << std::endl;
    las->PrintSolution(std::cout);
    */
    ///// update analysis with solution
    double * solution = NULL;
    las->GetSolution(solution);
    fem->UpdateDOFs(solution);
  }
}
