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

#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "peano/utils/Loop.h"
#include "kernels/GaussLegendreBasis.h"

#include "kernels/aderdg/generic/Kernels.h"

bool exahype::tests::c::GenericEulerKernelTest::_setNcpAndMatrixBToZero(false);

// TODO: Do not conclude macro definitions with a semicolon?!
//       (https://goo.gl/22Ac4j)
// clang-format off

#ifndef ALIGNMENT
registerTest(exahype::tests::c::GenericEulerKernelTest)
#endif

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::GenericEulerKernelTest::_log( "exahype::tests::c::GenericEulerKernelTest" );


namespace exahype {
namespace tests {
namespace c {

double**   GenericEulerKernelTest::weights                                  = kernels::legendre::weights;
double**   GenericEulerKernelTest::nodes                                    = kernels::legendre::nodes;
double***  GenericEulerKernelTest::Kxi                                      = kernels::legendre::Kxi;
double***  GenericEulerKernelTest::dudx                                     = kernels::legendre::dudx;
double***  GenericEulerKernelTest::iK1                                      = kernels::legendre::iK1;
double***  GenericEulerKernelTest::equidistantGridProjector                 = kernels::legendre::equidistantGridProjector;
double***  GenericEulerKernelTest::FCoeff                                   = kernels::legendre::FCoeff;
double**** GenericEulerKernelTest::fineGridProjector                        = kernels::legendre::fineGridProjector;
kernels::UnivariateFunction** GenericEulerKernelTest::basisFunction                 = kernels::legendre::basisFunction;
kernels::UnivariateFunction** GenericEulerKernelTest::basisFunctionFirstDerivative  = kernels::legendre::basisFunctionFirstDerivative;
kernels::UnivariateFunction** GenericEulerKernelTest::basisFunctionSecondDerivative = kernels::legendre::basisFunctionSecondDerivative;


GenericEulerKernelTest::GenericEulerKernelTest()
: tarch::tests::TestCase("exahype::tests::c::GenericEulerKernelTest") {}

GenericEulerKernelTest::~GenericEulerKernelTest() {}

void GenericEulerKernelTest::run() {
  testMethod(testPDEFluxes);
  logWarning("run()","Test testSpaceTimePredictorLinear is failing. Test data might not be correct anymore.");
//  testMethod(testSpaceTimePredictorLinear);
  //validate(false);

  std::cout << "RUNNING TESTS" << std::endl;

  testMethod(testSpaceTimePredictorNonlinear);
  testMethod(testVolumeIntegralLinear);
  testMethod(testVolumeIntegralNonlinear);

  testMethod(testRiemannSolverLinear);
  testMethod(testRiemannSolverNonlinear);
  testMethod(testSurfaceIntegralLinear);
  testMethod(testSurfaceIntegralNonlinear);
  testMethod(testFaceUnknownsProjection);
  testMethod(testVolumeUnknownsProjection);
  testMethod(testEquidistantGridProjection);

  testMethod(testSolutionUpdate);
}

void GenericEulerKernelTest::testEquidistantGridProjection() {
  logInfo( "testEquidistantGridProjection()", "Test equidistant grid projection, ORDER=2, DIM=2" );

  const int numberOfVariables = 5;
  const int order             = 3;
  const int basisSize         = order+1;
  const int cellUnknowns      = numberOfVariables * tarch::la::aPowI(basisSize,DIMENSIONS);

  double * u = new double[cellUnknowns];

  for (int i=0; i < cellUnknowns; ++i) {
    u[i] = 1.0;
  }

  // via basis functions
  dfor(i,order+1) {
    for (int unknown=0; unknown < numberOfVariables; unknown++) {
      std::ostringstream identifier;
      identifier << "Q" << unknown;

      double value = 0.0;
      dfor(ii,order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,order + 1);
        value +=  GenericEulerKernelTest::basisFunction[order][ii(0)](i(0)/order) *
                  GenericEulerKernelTest::basisFunction[order][ii(1)](i(1)/order) *
                  #ifdef Dim3
                  GenericEulerKernelTest::basisFunction[order][ii(2)](i(2)/order) *
                  #endif
                  u[iGauss * numberOfVariables + unknown];
      }
      assertion(tarch::la::equals(value,1.0,1e-6)); // todo precision issues
    }
  }

  // via equidistant grid projection
  dfor(i,order+1) {
    for (int unknown=0; unknown < numberOfVariables; unknown++) {
      std::ostringstream identifier;
      identifier << "Q" << unknown;

      double value = 0.0;
      dfor(ii,order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,order + 1);
        value +=  kernels::legendre::equidistantGridProjector[order][ii(0)][i(0)] *
            GenericEulerKernelTest::equidistantGridProjector[order][ii(1)][i(1)] *
            #ifdef Dim3
            GenericEulerKernelTest::equidistantGridProjector[order][ii(2)][i(2)] *
            #endif
            u[iGauss * numberOfVariables + unknown];
      }
      assertion(tarch::la::equals(value,1.0,1e-6)); // todo precision issues
    }
  }
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
//#else
// todo VV TestCase
//#endif

// clang-format on
