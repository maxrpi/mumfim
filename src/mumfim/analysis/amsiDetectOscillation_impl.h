#ifndef AMSI_DETECT_OSCILLATION_IMPL_H_
#define AMSI_DETECT_OSCILLATION_IMPL_H_
#include <amsiDetectOscillation.h>
#include <amsiNonlinearAnalysis.h>
#include <amsiVerbosity.h>
namespace amsi {
  template <typename T>
  struct DetectOscillationIteration : public Iteration {
    DetectOscillationIteration(T cnvg, Iteration* itr, unsigned int maxItr)
        : cnvg(cnvg), itr(itr), maxItr(maxItr)
    {
    }
    virtual void iterate() 
    {
      // so we want to subtract 1 to get the actual current iteration number
      if (itr->iteration() - 1 >= maxItr) {
        AMSI_V2(std::cout << "Solution has surpassed iteration cap of "
                          << itr->iteration() - 1 << "\n";)
        fail = true;
        Iteration::iterate();
      }
    }
    protected:
    T cnvg;
    Iteration* itr;
    unsigned int maxItr;
  };
  template <typename T>
  struct DetectOscillationPrevNorm : public Iteration {
    DetectOscillationPrevNorm(T cnvg, Iteration* itr, double prevNormFactor)
        : cnvg(cnvg), itr(itr), prevNormFactor(prevNormFactor)
    {
    }
    virtual void iterate() 
    {
      if (cnvg->getPrevNorm() * prevNormFactor < cnvg->getCurrNorm()) {
        AMSI_V2(std::cout << "Previous Norm *" << prevNormFactor
                          << " > Current Norm\n";
                std::cout << "Previous Norm: "
                          << cnvg->getPrevNorm() * prevNormFactor << " "
                          << "Current Norm: " << cnvg->getCurrNorm() << "\n";)
        fail = true;
      }
    }
    protected:
    T cnvg;
    Iteration* itr;
    double prevNormFactor;
  };
  template <typename T>
  Iteration* createOscillationDetection(DetectOscillationType type, T cnvg,
                                        Iteration* itr,
                                        unsigned int maxIteration,
                                        double prevNormFactor)
  {
    switch (type) {
      case DetectOscillationType::IterationOnly:
        return new DetectOscillationIteration<T>(cnvg, itr, maxIteration);
      case DetectOscillationType::PrevNorm:
        return new DetectOscillationPrevNorm<T>(cnvg, itr, prevNormFactor);
      case DetectOscillationType::IterationPrevNorm:
        Iteration* dtct_itr =
            new DetectOscillationIteration<T>(cnvg, itr, maxIteration);
        Iteration* dtct_prev_nrm =
            new DetectOscillationPrevNorm<T>(cnvg, itr, prevNormFactor);
        std::vector<Iteration*> itr_stps;
        itr_stps.push_back(dtct_itr);
        itr_stps.push_back(dtct_prev_nrm);
        return new MultiIteration(itr_stps.begin(), itr_stps.end());
    }
    return NULL;
  }
}  // namespace amsi
#endif
