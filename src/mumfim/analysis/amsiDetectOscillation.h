#ifndef AMSI_DETECT_OSCILLATION_H_
#define AMSI_DETECT_OSCILLATION_H_
#include <amsiNonlinearAnalysis.h>
namespace amsi {
  /* Method for detecting oscillations */
  enum class DetectOscillationType {
    /* Detects if the solution has surpassed the maximum number of oscillations
     */
    IterationOnly,
    /* Detects if the previous norm is bigger than the current norm times some
       factor*/
    PrevNorm,
    /* Checks that the iteration is below the max number of iterations and the
       current norm times some factor is less than the previous norm */
    IterationPrevNorm
  };
  // the convergence type must implement getPrevNorm and getCurrNorm functions
  template <typename T>
  Iteration* createOscillationDetection(DetectOscillationType type, T cnvg,
                                        amsi::Iteration* itr,
                                        unsigned int maxIteration,
                                        double prevNormFactor);
}  // namespace amsi
#include "amsiDetectOscillation_impl.h"
#endif
