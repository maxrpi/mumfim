#ifndef AMSI_NONLINEAR_ANALYSIS_H_
#define AMSI_NONLINEAR_ANALYSIS_H_
#include <amsiVerbosity.h>
#include <cassert>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
namespace amsi {
  class Iteration;
  class Convergence;
  /**
   * Perform a numerical solve using the supplied iteration
   *  and convergence objects. The iteration object
   *  constitutes a single iteration in a nonlinear solve.
   *  The convergence object informs whether the underlying
   *  simulation has converged to a solution for the current
   *  nonlinear solve.
   *  If the simulation is linear, simply pass a convergence
   *  object that returns true and an iteration object that
   *  performs the linear solve.
   */
  bool numericalSolve(Iteration* it, Convergence* cn);
  struct PerIter {
    virtual void iter() = 0;
  };
  struct PerStep {
    virtual void step() = 0;
  };
  /**
   * An iteration object that represents a single iteration
   *  for a numerical simulation. The base class merely tracks
   *  the iteration count as a single linear solve takes place.
   *  Classes which inherit from the Iteration should form the
   *  linear system for a simulation and perform a solve of
   *  of the linear system. A suitable convergence operator
   *  to detect whether the state of the simulation after
   *  the linear system constitutes a solved state.
   */
  class Iteration {
    protected:
    unsigned long itr;
    bool fail;

    public:
    Iteration() : itr(0), fail(false) {}
    virtual ~Iteration() {}
    virtual void iterate() { ++itr; }
    unsigned long iteration() const { return itr; }
    virtual bool failed() { return fail; }
    virtual void reset() { itr = 0; }
  };
  /*
   * A class that allows the construction of an interation out of list of
   * iteration objects. This is used in the same way as MultiConverence. It is
   * also similar to ModularIteration, however since it contains full iteration
   * objects it gives more control over the ordering of things.
   * MultiIteration does not take ownership of any sub iterations and expects
   * the user to delete them after use
   *
   */
  class MultiIteration : public Iteration {
    protected:
    std::vector<Iteration*> itrs;
    public:
    template <typename I>
    MultiIteration(I bgn, I end) : Iteration(), itrs()
    {
      std::copy(bgn, end, std::back_inserter(itrs));
    }
    // do not assume ownership of sub iterations e.g. do not delete them
    // we do this for 2 reasons.
    // 1. we may want to access the subiteration after the multiiteration is
    // complete (unlikely)
    // 2. we do not want to assume all sub iterations are allocated on the heap
    virtual ~MultiIteration() {}
    virtual void iterate()
    {
      Iteration::iterate();
      for (auto itr = itrs.begin(); itr != itrs.end(); ++itr) {
        (*itr)->iterate();
      }
    }
    virtual bool failed()
    {
      for (auto itr = itrs.begin(); itr != itrs.end(); ++itr) {
        if ((*itr)->failed()) {
          fail = true;
          return true;
        }
      }
      return false;
    }
    virtual void reset()
    {
      for (auto itr = itrs.begin(); itr != itrs.end(); ++itr) {
        (*itr)->reset();
      }
      // reset the global iteration too
      this->itr = 0;
    }
    virtual void addIteration(Iteration* itr) { itrs.push_back(itr); }
  };
  /**
   * An operation which determines whether the
   *  modeled operation has converged.
   */
  class Convergence {
    public:
    Convergence() : fail(false) {}
    virtual bool converged() = 0;
    virtual bool failed() { return fail; }
    virtual ~Convergence(){};

    protected:
    bool fail;
  };
  /**
   * The updating convergence class allows the internal
   *  epsilon value to change depending on the state of
   *  the simulation.
   * @tparam V the current convergence value object,
   *            has a double operator()() function which
   *            gives the current convergence value.
   * @tparam E the current epsilon value object,
   *            has a double operator()() function which
   *            gives the current epsilon value.
   * @tparam R the current reference value object,
   *            has a double operator()() function which
   *            gives the current reference value.
   */
  template <typename V, typename E, typename R>
  class UpdatingConvergence : public Convergence {
    protected:
    Iteration* itr;
    double cvg_vl;
    double prev_vl;
    double eps;
    double ref_vl;
    V cvg_gen;
    E eps_gen;
    R ref_gen;
    bool failOnNonConvergence;

    public:
    // prev_vl is initialized greater than cvg_vl so oscillation detection have
    // problems
    UpdatingConvergence(Iteration* it, V val_gen, E eps_gen, R ref_gen,
                        bool failOnNonConvergence = false)
        : itr(it)
        , cvg_vl(std::numeric_limits<double>::max())
        , prev_vl(std::numeric_limits<double>::max())
        , eps(1e-16)
        , ref_vl(std::numeric_limits<double>::max())
        , cvg_gen(val_gen)
        , eps_gen(eps_gen)
        , ref_gen(ref_gen)
        , failOnNonConvergence(failOnNonConvergence)
    {
    }
    ~UpdatingConvergence() {}
    /**
     * Update the internal values used to determine
     *  convergence, the test value, epsilon, and
     *  reference are updated from their respective
     *  objects.
     */
    virtual void update()
    {
      prev_vl = cvg_vl;
      cvg_vl = (*cvg_gen)();
      eps = (*eps_gen)(itr->iteration());
      ref_vl = (*ref_gen)();
    }
    double getPrevNorm() const { return prev_vl; }
    double getCurrNorm() const { return cvg_vl; }
    /**
     * Determine whether the modeled operation has
     *  converged.
     * @return \f$ test < \eps * ref \f$
     */
    virtual bool converged()
    {
      update();
      bool cvrgd = cvg_vl <= eps * ref_vl;
      AMSI_V1(std::cout << "convergence criteria: \n"
                        << "\t" << cvg_vl << " < " << eps << " * " << ref_vl
                        << "\n"
                        << "\t" << cvg_vl << " < " << eps * ref_vl << "\n"
                        << "\t" << (cvrgd ? "TRUE" : "FALSE") << std::endl;)
      if (!cvrgd && failOnNonConvergence) {
        this->fail = true;
      }
      return cvrgd;
    }
  };
  /**
   * A convergence class that wraps multiple convergence
   *  objects. Useful for composing simple convergence
   *  operations in lieu of implementing a bespoke
   *  combined class. Cannot modify the convergence
   *  operations after construction.
   * This is kind of a variation on the Composite pattern.
   * MultiConvergence does not take ownership of any sub convergence and expects
   * the user to delete them after use
   */
  class MultiConvergence : public Convergence {
    private:
    std::vector<Convergence*> cvgs;

    public:
    template <typename I>
    MultiConvergence(I bgn, I end) : Convergence(), cvgs()
    {
      std::copy(bgn, end, std::back_inserter(cvgs));
    }
    // do not assume ownership of sub iterations e.g. do not delete them
    // we do this for 2 reasons.
    // 1. we may want to access the subiteration after the multiiteration is
    // complete (unlikely)
    // 2. we do not want to assume all sub iterations are allocated on the heap
    virtual ~MultiConvergence() {}
    virtual bool converged()
    {
      for (auto cvg = cvgs.begin(); cvg != cvgs.end(); ++cvg)
        if (!(*cvg)->converged()) return false;
      return true;
    }
    virtual bool failed()
    {
      for (auto cvg = cvgs.begin(); cvg != cvgs.end(); ++cvg)
        if ((*cvg)->failed()) {
          fail = true;
          return true;
        }
      return false;
    }
  };
  /**
   * A convergence class to use for linear
   *  simulations, since no iteration is required
   *  for linear simulations, this simply returns
   *  true.
   */
  class LinearConvergence : public Convergence {
    public:
    virtual bool converged() { return true; }
  };
  extern LinearConvergence
      linear_convergence;  // should be const but converged() isn't const
  // an iteration which will stop when the max number of iterations is reached
  struct StopAtMaxIters : public Iteration {
    StopAtMaxIters(unsigned int maxItr) : maxItr(maxItr) {}
    virtual void iterate()
    {
      if ((maxItr == 0) || this->iteration() >= maxItr - 1) {
        AMSI_V2(std::cout << "Solution has surpassed iteration cap of "
                          << this->iteration() << "\n";)
        fail = true;
      }
      Iteration::iterate();
    }

    protected:
    unsigned int maxItr;
  };
  // convenience function for type inference with the updating convergencce
  template <typename V, typename E, typename R>
  std::unique_ptr<UpdatingConvergence<V, E, R> > createUpdatingConvergence(
      Iteration* it, V val_gen, E eps_gen, R ref_gen,
      bool failOnConvergence = false)
  {
    return std::make_unique<UpdatingConvergence<V,E,R>>(it, val_gen, eps_gen, ref_gen);
  }
}  // namespace amsi
#endif
