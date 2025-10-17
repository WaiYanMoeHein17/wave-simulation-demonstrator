/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "EulerSolver_ADERDG.h"
#include "EulerSolver_ADERDG_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <algorithm>

#include <string>

#include <math.h>

#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h"

#include "kernels/GaussLegendreBasis.h"

tarch::logging::Log Euler::EulerSolver_ADERDG::_log("Euler::EulerSolver_ADERDG");

Euler::EulerSolver_ADERDG::Reference Euler::EulerSolver_ADERDG::ReferenceChoice = Euler::EulerSolver_ADERDG::Reference::EntropyWave;

void Euler::EulerSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  if (constants.isValueValidString("reference")) {
    std::string reference = constants.getValueAsString("reference");
    
    if (reference.compare("entropywave")==0) {
      ReferenceChoice = Reference::EntropyWave;
    }
    else if (reference.compare("rarefactionwave")==0) {
      ReferenceChoice = Reference::RarefactionWave;
    }
    else if (reference.compare("sod")==0){
      ReferenceChoice = Reference::SodShockTube;
    }
    else if (reference.compare("explosion")==0){
      ReferenceChoice = Reference::SphericalExplosion;
    }
    else {
      logError("init(...)","do not recognise value '"<<reference<<"' for constant 'reference'. Use either 'entropywave', "
              "'rarefactionwave', 'sod', or 'explosion'.");
      std::abort();
    }
    logInfo("init(...)","use initial condition '" << reference << "'.");
  } else {
    logInfo("init(...)","use initial condition 'entropyWave' (default value).");
  }
}

void Euler::EulerSolver_ADERDG::flux(const double* const Q, double** const F) {
  #ifdef SymbolicVariables
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double gamma = 1.4;
  const double irho = 1./vars.rho();
  const double p = (gamma-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
  #else // SymbolicVariables
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  // col 1
  F[0][0] = Q[1];
  F[0][1] = irho*Q[1]*Q[1] + p;
  F[0][2] = irho*Q[2]*Q[1];
  F[0][3] = irho*Q[3]*Q[1];
  F[0][4] = irho*(Q[4]+p)*Q[1];

  // col 2
  F[1][0] = Q[2];
  F[1][1] = irho*Q[1]*Q[2];
  F[1][2] = irho*Q[2]*Q[2] + p;
  F[1][3] = irho*Q[3]*Q[2];
  F[1][4] = irho*(Q[4]+p)*Q[2];

  #if DIMENSIONS==3
  // col 3
  F[2][0] = Q[3];
  F[2][1] = irho*Q[1]*Q[3];
  F[2][2] = irho*Q[2]*Q[3];
  F[2][3] = irho*Q[3]*Q[3] + p;
  F[2][4] = irho*(Q[4]+p)*Q[3];
  #endif
  #endif
}

void Euler::EulerSolver_ADERDG::eigenvalues(const double* const Q,
    const int direction,
    double* const lambda) {
  #ifdef SymbolicVariables
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double gamma = 1.4;
  const double irho = 1./vars.rho();
  const double p = (gamma-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[direction + 1] * irho;
  double c   = std::sqrt(gamma * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
  #else // SymbolicVariables
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  const double u_n = Q[direction + 1] * irho;
  const double c   = std::sqrt(gamma * p * irho);

  std::fill_n(lambda,5,u_n);
  lambda[0] -= c;
  lambda[4] += c;
  #endif
}

void Euler::EulerSolver_ADERDG::entropyWave(const double* const x,double t, double* const Q) {
  const double gamma         = 1.4;
  constexpr double width     = 0.3;
  constexpr double amplitude = 0.3;

  #if DIMENSIONS==2
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5);
  #else
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1],x[2]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.0,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5,0.5);
  #endif
  const double distance  = tarch::la::norm2( xVec - x0 - v0 * t );

  Q[0] = 0.5 + amplitude * std::exp(-distance / std::pow(width, DIMENSIONS));
  Q[1] = Q[0] * v0[0];
  Q[2] = Q[0] * v0[1];
  Q[3] = 0.0;
  // total energy = internal energy + kinetic energy
  const double p = 1.;
  Q[4] = p / (gamma-1)  +  0.5*Q[0] * (v0[0]*v0[0]+v0[1]*v0[1]); // v*v; assumes: v0[2]=0
}

void Euler::EulerSolver_ADERDG::sodShockTube(const double* const x, const double t, double* const Q) {
  // Initial data
  constexpr double gamma     =1.39999999999999991118;
  constexpr double x_0       =0.50000000000000000000;

  constexpr double rho_5     =0.12500000000000000000; // right states
  constexpr double P_5       =0.10000000000000000555;
  constexpr double u_5       =0.00000000000000000000;
  constexpr double rho_1     =1.00000000000000000000; // left states
  constexpr double P_1       =1.00000000000000000000;
  constexpr double u_1       =0.00000000000000000000;

  // Sound speed
  constexpr double cs_1       =1.18321595661992318149;

  // Contact left
  constexpr double rho_3     =0.42631942817849538541;
  constexpr double P_3       =0.30313017805064701449;
  constexpr double u_3       =0.92745262004895057117;
  constexpr double cs_3      =0.99772543261013335592;

  // Contact right
  constexpr double rho_4     =0.26557371170530713611;
  constexpr double P_4       =0.30313017805064701449;
  constexpr double u_4       =0.92745262004895057117;

  // Shock
  constexpr double u_shock   =1.75215573203017838111;

  // Key Positions
  const double x_4 = x_0 + u_shock * t;      // position of shock
  const double x_3 = x_0 + u_3 * t;          // position of contact discontinuity
  const double x_2 = x_0 + (u_3 - cs_3) * t; // foot of rarefaction wave
  const double x_1 = x_0 - cs_1 * t;         // head of rarefaction wave

  double p = 0; // pressure
  Q[2] = 0; // y velocity
  Q[3] = 0; // z velocity
  if (tarch::la::equals(t,0.0)) {
    if (x[0] < x_0) {
      Q[0] = rho_1;
      Q[1] = Q[0] * u_1;
      p    = P_1;
    } else {
      Q[0] = rho_5;
      Q[1] = Q[0] * u_5;
      p    = P_5;
    }
    Q[4] = p/(gamma-1) + 0.5 / Q[0] * (Q[1]*Q[1]); // j*j, j=rho*v !!! ; assumes: Q[1+i]=0, i=1,2.

  } else {
    if (x[0] < x_1) {
      Q[0] = rho_1;
      Q[1] = Q[0] * u_1;
      p    = P_1;
    } else if (x_1 <= x[0] && x[0] < x_2) {
      // rarefaction wave
      const double u      = 2.0 / (gamma+1) * (cs_1 + (x[0] - x_0) / t);
      const double factor = 1.0 - 0.5*(gamma-1)*u / cs_1;
      Q[0] = rho_1 * std::pow( factor, 2/(gamma-1) );
      Q[1] = Q[0]  * u;
      p    = P_1   * std::pow( factor, 2.0*gamma/(gamma-1) );
    } else if (x_2 <= x[0] && x[0] < x_3) {
      Q[0] = rho_3;
      Q[1] = Q[0] * u_3;
      p    = P_3;
    } else if (x_3 <= x[0] && x[0] < x_4) {
      Q[0] = rho_4;
      Q[1] = Q[0] * u_4;
      p    = P_4;
    } else if (x_4 <= x[0]) {
      Q[0] = rho_5;
      Q[1] = Q[0] * u_5;
      p    = P_5;
    }
    // total energy = internal energy + kinetic energy
    Q[4] = p/(gamma-1) + 0.5 / Q[0] * (Q[1]*Q[1]); // j*j, j=rho*v !!! ; assumes: Q[1+i]=0, i=1,2.
  }
}

void Euler::EulerSolver_ADERDG::sphericalExplosion(const double* const x,double t, double* const Q) {
  constexpr double x0[3]   = {0.5, 0.5, 0.5};
  constexpr double radius  = 0.25;
  constexpr double radius2 = radius*radius;

  // Velocities are set to zero (initially).
  if (tarch::la::equals(t,0.0)) {
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = 0.0;
    #if DIMENSIONS==2
    // Circular shaped pressure jump at centre of domain.
    if((x[0]-x0[0])*(x[0]-x0[0]) + (x[1]-x0[1])*(x[1]-x0[1]) < radius2) {
      Q[0] = 1.0;
      Q[4] = 1.0;
    } else {
      Q[0] = 0.125;
      Q[4] = 0.1; // total energy
    }
    #else
    // Circular shaped pressure jump at centre of domain.
    if((x[0]-x0[0])*(x[0]-x0[0]) + (x[1]-x0[1])*(x[1]-x0[1]) + (x[2]-x0[2])*(x[2]-x0[2]) < radius2) {
      Q[0] = 1.0;
      Q[4] = 1.0;
    } else {
      Q[0] = 0.125;
      Q[4] = 0.1; // total energy
    }
    #endif
  } else {
    std::fill_n(Q, NumberOfVariables, 0.0);
    // We then compute the norm in our error writers for t>0.
  }
}

void Euler::EulerSolver_ADERDG::rarefactionWave(const double* const x,double t, double* const Q) {
  constexpr double gamma = 1.4;
  constexpr double width = 0.25;
  constexpr double x0[3] = { 0.5, 0.5, 0.5 };

  if (tarch::la::equals(t,0.0)) {
    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    #if DIMENSIONS==2
    const double norm2Squared = (x[0]-x0[0])*(x[0]-x0[0]) + (x[1]-x0[1])*(x[1]-x0[1]);
    #else
    const double norm2Squared = (x[0]-x0[0])*(x[0]-x0[0]) + (x[1]-x0[1])*(x[1]-x0[1]) + (x[2]-x0[2])*(x[2]-x0[2]);
    #endif
    Q[4] = 1. / (gamma - 1) + // pressure is set to one
        exp(-std::sqrt(norm2Squared) / pow(width, DIMENSIONS)) * 2;
  }
  else {
    std::fill_n(Q, NumberOfVariables, 0.0);
    // We then compute the norm in our error writers for t>0.
  }
}

void Euler::EulerSolver_ADERDG::referenceSolution(const double* const x,double t, double* const Q) {
  switch (ReferenceChoice) {
  case Reference::SodShockTube:
    sodShockTube(x,t,Q);
    break;
  case Reference::EntropyWave:
    entropyWave(x,t,Q);
    break;
  case Reference::SphericalExplosion:
    sphericalExplosion(x,t,Q);
    break;
  case Reference::RarefactionWave:
    rarefactionWave(x,t,Q);
    break;
  }
}

void Euler::EulerSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t, 0.0)) {
    referenceSolution(x,0.0,Q);
  }
  //referenceSolution(x,0.0,Q);

}

exahype::solvers::Solver::RefinementControl
Euler::EulerSolver_ADERDG::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  double maxE = -std::numeric_limits<double>::max();
  double minE = +std::numeric_limits<double>::max();
  
  const int nodes = std::pow((Order+1), DIMENSIONS);
  for(int i=0;i<nodes;i++) {
    maxE = std::max(maxE, luh[i*NumberOfVariables+NumberOfVariables-1]);
    minE = std::min(minE, luh[i*NumberOfVariables+NumberOfVariables-1]);
  }

  if ( maxE/minE > 1.05 ) {
    return RefinementControl::Refine;
  } else if ( level > getCoarsestMeshLevel() ) {
    return RefinementControl::Erase;
  } else {
    return RefinementControl::Keep;
  }
}

void Euler::EulerSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  switch (ReferenceChoice) {
  case Reference::SphericalExplosion:
  case Reference::RarefactionWave:
  case Reference::SodShockTube: { // wall boundary conditions
    std::copy_n(stateIn, NumberOfVariables, stateOut);
    stateOut[1+direction] =  -stateOut[1+direction];
    double _F[3][NumberOfVariables]={0.0};
    double* F[3] = {_F[0], _F[1], _F[2]};
    flux(stateOut,F);
    std::copy_n(F[direction], NumberOfVariables, fluxOut);
  } break;
  case Reference::EntropyWave: {// Dirichlet conditions
    double Q[NumberOfVariables]     = {0.0};
    double _F[3][NumberOfVariables] = {0.0};
    double* F[3] = {_F[0],_F[1],_F[2]};

    // initialise
    std::fill_n(stateOut, NumberOfVariables, 0.0);
    std::fill_n(fluxOut,  NumberOfVariables, 0.0);
    for (int i=0; i<Order+1; i++) {
      const double ti = t + dt * kernels::legendre::nodes[Order][i];
      referenceSolution(x,ti,Q);
      flux(Q,F);
      for (int v=0; v<NumberOfVariables; v++) {
        stateOut[v] += Q[v]            * kernels::legendre::weights[Order][i];
        fluxOut[v]  += F[direction][v] * kernels::legendre::weights[Order][i];
      }
    }
  } break;
  }
}

void Euler::EulerSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
    std::copy_n(Q, NumberOfVariables, observables);
}


bool Euler::EulerSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  if (observablesMin[0] <= 0.0) return false;
  if (observablesMin[4] < 0.0) return false;
  return true;
}
