#include "ErrorPlotter.h"

#include "EulerSolver_ADERDG.h"

Euler::ErrorPlotter::ErrorPlotter(EulerSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}


Euler::ErrorPlotter::~ErrorPlotter() {
  // @todo Please insert your code here
}


void Euler::ErrorPlotter::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::ErrorPlotter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::ErrorPlotter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  constexpr int numberOfVariables = AbstractEulerSolver_ADERDG::NumberOfVariables;

  double QAna[numberOfVariables];
  EulerSolver_ADERDG::referenceSolution(x.data(),timeStamp,QAna);

  for (int i=0; i<numberOfVariables; i++){ 
    outputQuantities[i] =std::abs( QAna[i] - Q[i] );
  }
  for (int i=0; i<numberOfVariables; i++){ 
    outputQuantities[i+numberOfVariables] = Q[i];
  } 
  for (int i=0; i<numberOfVariables; i++){ 
    outputQuantities[i+2*numberOfVariables] = QAna[i];
  }
}


