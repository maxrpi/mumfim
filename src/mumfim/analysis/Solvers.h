#ifndef FEASOLVERS_H_
#define FEASOLVERS_H_
#include "amsiFEA.h"
#include "amsiLAS.h"
#include "amsiNonlinearAnalysis.h"
#include <string>
namespace amsi
{
  Iteration * buildLinearFEMIteration(FEAStep * f, LAS * l);
  void LinearSolver(FEAStep *,LAS*);
  //void NewtonSolver(FEAStep*,LAS*,int,double,double,double&);
}
#endif
