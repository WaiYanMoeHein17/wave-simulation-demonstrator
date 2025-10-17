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

#include "exahype/tests/kernels/c/ElasticityKernelTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "kernels/GaussLegendreBasis.h"
#include "peano/utils/Loop.h"

#include "kernels/aderdg/generic/Kernels.h"

// TODO: Do not conclude macro definitions with a semicolon?!
//       (https://goo.gl/22Ac4j)
// clang-format off
#ifndef ALIGNMENT
registerTest(exahype::tests::c::ElasticityKernelTest)
#endif

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::ElasticityKernelTest::_log( "exahype::tests::c::ElasticityKernelTest" );


namespace exahype {
namespace tests {
namespace c {

double**   ElasticityKernelTest::weights                                  = kernels::lobatto::weights;
double**   ElasticityKernelTest::nodes                                    = kernels::lobatto::nodes;
double***  ElasticityKernelTest::Kxi                                      = kernels::lobatto::Kxi;
double***  ElasticityKernelTest::dudx                                     = kernels::lobatto::dudx;
double***  ElasticityKernelTest::iK1                                      = kernels::lobatto::iK1;
double***  ElasticityKernelTest::equidistantGridProjector                 = kernels::lobatto::equidistantGridProjector;
double***  ElasticityKernelTest::FCoeff                                   = kernels::lobatto::FCoeff;
double**** ElasticityKernelTest::fineGridProjector                        = kernels::lobatto::fineGridProjector;
kernels::UnivariateFunction** ElasticityKernelTest::basisFunction                 = kernels::lobatto::basisFunction;
kernels::UnivariateFunction** ElasticityKernelTest::basisFunctionFirstDerivative  = kernels::lobatto::basisFunctionFirstDerivative;
kernels::UnivariateFunction** ElasticityKernelTest::basisFunctionSecondDerivative = kernels::lobatto::basisFunctionSecondDerivative;


ElasticityKernelTest::ElasticityKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::ElasticityKernelTest") {}

ElasticityKernelTest::~ElasticityKernelTest() {}

void ElasticityKernelTest::run() {
  _log.info("ElasticityKernelTest::run()", "ElasticityKernelTest is active");
  // TODO(Dominic): Assess
//  testMethod(testPDEFluxes);
  logWarning("run()","Test testSpaceTimePredictorLinear is disabled! Test data might be outdated.");
//  testMethod(testSpaceTimePredictorLinear); // TODO(Dominic): Fix
//  testMethod(testSpaceTimePredictorNonlinear); // OPTIONAL
  testMethod(testVolumeIntegralLinear);
//  testMethod(testVolumeIntegralNonlinear); // OPTIONAL
  testMethod(testSurfaceIntegralLinear);
//  testMethod(testSurfaceIntegralNonlinear) // OPTIONAL
  logWarning("run()","Test testRiemannSolverLinear is disabled!");
  testMethod(testRiemannSolverLinear); // TODO(Dominic): Fix
//  testMethod(testRiemannSolverNonlinear);// OPTIONAL
//  testMethod(testFaceUnknownsProjection); // OPTIONAL
//  testMethod(testVolumeUnknownsProjection); // OPTIONAL
//  testMethod(testEquidistantGridProjection); // OPTIONAL
//
//  testMethod(testSolutionUpdate);
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
