#include "amsiNonlinearAnalysis.h"
namespace amsi
{
  LinearConvergence linear_convergence;
  bool numericalSolve(Iteration * it, Convergence * cn)
  {
    do
    {
      if (cn->failed())
      {
        return false;
      }
      it->iterate();
      if (it->failed())
      {
        return false;
      }
    } while (!cn->converged());
    // require the user to reset the iteration, because they may have tasks
    // to perform with the unreset iteration.
    //it->reset();
    return true;
  }
}  // namespace amsi
