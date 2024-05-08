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
  class LinearIteration : public Iteration
  {
  private:
  FEAStep * fem;
    LAS * las;
  public:
    LinearIteration(FEAStep * f, LAS * l)
      : fem(f)
      , las(l)
    { }
    void iterate()
    {
      fem->ApplyBC_Dirichlet();
      fem->RenumberDOFs();
      int gbl, lcl, off;
      fem->GetDOFInfo(gbl,lcl,off);
      las->Reinitialize(lcl,gbl,off);
      las->Zero();
      fem->Assemble(las);
      las->solve();
      double * s = NULL;
      las->GetSolution(s);
      fem->UpdateDOFs(s);
      Iteration::iterate();
    }
  };
  Iteration * buildLinearFEMIteration(FEAStep * f, LAS * l)
  {
    return new LinearIteration(f,l);
  }
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
  void NewtonSolver(FEAStep * fem,
                    LAS * las,
                    int iteration_cap,
                    double epsilon,
                    double,
                    double & residual_norm)
  {
    int global_dof_count, local_dof_count, first_local_dof;
    int newton_iteration = 0;
    while(true)
    {
      AMSI_V1(std::cout << "Newton iteration " << newton_iteration << ":" << std::endl;)
      fem->ApplyBC_Dirichlet();
      fem->RenumberDOFs();
      fem->GetDOFInfo(global_dof_count,local_dof_count,first_local_dof);
      las->Reinitialize(local_dof_count,global_dof_count,first_local_dof);
      las->Zero();
      fem->Assemble(las);
      double nrm = 0.0;
      las->GetVectorNorm(nrm);
      las->solve();
      double * solution = NULL;
      las->GetSolution(solution);
      fem->UpdateDOFs(solution);
      las->GetVectorNorm(residual_norm);
      if (residual_norm < epsilon)
        break;
      if (newton_iteration > iteration_cap)
        break;
      fem->Adapt();
      newton_iteration++;
    }
  }
}
