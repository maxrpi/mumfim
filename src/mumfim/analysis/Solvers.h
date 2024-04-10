#ifndef FEASOLVERS_H_
#define FEASOLVERS_H_
#include "amsiFEA.h"
#include "amsiLAS.h"
#include "amsiNonlinearAnalysis.h"
#include <string>
namespace amsi
{
  Iteration * buildLinearFEMIteration(FEA * f, LAS * l);
  void LinearSolver(FEA*,LAS*);
  //void NewtonSolver(FEA*,LAS*,int,double,double,double&);
}
#endif
