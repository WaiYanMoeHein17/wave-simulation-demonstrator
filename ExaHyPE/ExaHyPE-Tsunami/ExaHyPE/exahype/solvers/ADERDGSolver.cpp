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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Guera, Leonhard Rannabauer
 **/
#include "exahype/solvers/ADERDGSolver.h"

#include <limits>
#include <iomanip>
#include <chrono>
#include <algorithm> // copy_n

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "tarch/la/VectorVectorOperations.h"
#include "tarch/multicore/Lock.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/heap/CompressedFloatingPointNumbers.h"
#include "peano/datatraversal/TaskSet.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "kernels/KernelUtils.h"

#include "tarch/multicore/Jobs.h"
#include "tarch/la/Vector.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

#ifdef USE_ITAC
#include "VT.h"
#endif

namespace {
  constexpr const char* tags[]{"solutionUpdate",
                             "volumeIntegral",
                             "surfaceIntegral",
                             "riemannSolver",
                             "spaceTimePredictor",
                             "stableTimeStepSize",
                             "solutionAdjustment",
                             "faceUnknownsProlongation",
                             "faceUnknownsRestriction",
                             "volumeUnknownsProlongation",
                             "volumeUnknownsRestriction",
                             "boundaryConditions",
                             "deltaDistribution"
                             };
}

#ifdef USE_ITAC
int exahype::solvers::ADERDGSolver::adjustSolutionHandle                  = 0;
int exahype::solvers::ADERDGSolver::fusedTimeStepBodyHandle               = 0;
int exahype::solvers::ADERDGSolver::fusedTimeStepBodyHandleSkeleton       = 0;
int exahype::solvers::ADERDGSolver::predictorBodyHandle                   = 0;
int exahype::solvers::ADERDGSolver::predictorBodyHandleSkeleton           = 0;
int exahype::solvers::ADERDGSolver::updateBodyHandle                      = 0;
int exahype::solvers::ADERDGSolver::mergeNeighboursHandle                 = 0;
int exahype::solvers::ADERDGSolver::prolongateFaceDataToVirtualCellHandle = 0;
int exahype::solvers::ADERDGSolver::restrictToTopMostParentHandle         = 0;
#endif

tarch::logging::Log exahype::solvers::ADERDGSolver::_log( "exahype::solvers::ADERDGSolver");

// communication status
int exahype::solvers::ADERDGSolver::LeafCommunicationStatus                             = 2;
int exahype::solvers::ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication = 1;
// augmentation status
// On-the fly erasing seems to work with those values
int exahype::solvers::ADERDGSolver::MaximumAugmentationStatus                   = 3;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForVirtualRefining = 2;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForRefining        = 1;

/**
 * static constexpr need to declared again when following a
 * C++ standard before C++17.
 */
constexpr int exahype::solvers::ADERDGSolver::EmptyStatus;
constexpr int exahype::solvers::ADERDGSolver::BoundaryStatus;
constexpr int exahype::solvers::ADERDGSolver::Pending;
constexpr int exahype::solvers::ADERDGSolver::Erase; 
constexpr int exahype::solvers::ADERDGSolver::Keep;

tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::RestrictionSemaphore;

tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::CoarseGridSemaphore;

int exahype::solvers::ADERDGSolver::computeWeight(const int cellDescriptionsIndex) {
  if ( ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    int result = 0;
    for ( CellDescription& cellDescription : getCellDescriptions(cellDescriptionsIndex) ) {
      result += ( cellDescription.getType()==CellDescription::Type::Leaf ) ?  1 : 0;
    }
    return result;
  }
  else return 0;
}

/**
 * Returns the ADERDGCellDescription heap vector
 * at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::Heap::HeapEntries& exahype::solvers::ADERDGSolver::getCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  return Heap::getInstance().getData(cellDescriptionsIndex);
}

/**
 * Returns the ADERDGCellDescription with index \p element
 * in the heap vector at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::CellDescription& exahype::solvers::ADERDGSolver::getCellDescription(
    const int cellDescriptionsIndex,
    const int element) {
  assertion2(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex,element);
  assertion2(element>=0,cellDescriptionsIndex,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),cellDescriptionsIndex,element);

  return Heap::getInstance().getData(cellDescriptionsIndex)[element];
}

bool exahype::solvers::ADERDGSolver::communicateWithNeighbour(const CellDescription& cellDescription,const int faceIndex) {
  assertion1(cellDescription.getType()!=CellDescription::Type::Leaf ||
            cellDescription.getCommunicationStatus()==LeafCommunicationStatus,cellDescription.toString());
  return
      (cellDescription.getCommunicationStatus()                  == LeafCommunicationStatus &&
      cellDescription.getFacewiseCommunicationStatus(faceIndex)  >= MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getFacewiseAugmentationStatus(faceIndex)   <  MaximumAugmentationStatus)
      ||
      (cellDescription.getFacewiseCommunicationStatus(faceIndex) == LeafCommunicationStatus &&
      cellDescription.getCommunicationStatus()                   >= MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getAugmentationStatus()                    <  MaximumAugmentationStatus);
}

void exahype::solvers::ADERDGSolver::prefetchFaceData(CellDescription& cellDescription,const int faceIndex) {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  auto* solver = solvers::RegisteredSolvers[cellDescription.getSolverNumber()];
  int dataPerFace = 0;
  int dofPerFace  = 0;
  switch (solver->getType()) {
    case Solver::Type::ADERDG:
      dataPerFace = static_cast<ADERDGSolver*>(solver)->getBndFaceSize();
      dofPerFace  = static_cast<ADERDGSolver*>(solver)->getBndFluxSize();
      break;
    case Solver::Type::LimitingADERDG:
      dataPerFace = static_cast<LimitingADERDGSolver*>(solver)->getSolver()->getBndFaceSize();
      dofPerFace  = static_cast<LimitingADERDGSolver*>(solver)->getSolver()->getBndFluxSize();
      break;
    default:
      break;
  }

  double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + faceIndex * dataPerFace;
  double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation())           + faceIndex * dofPerFace;

  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}

exahype::solvers::ADERDGSolver::ADERDGSolver(
    const std::string& identifier,
    const int numberOfVariables,
    const int numberOfParameters,
    const int numberOfGlobalObservables,
    const int basisSize,
    const double maximumMeshSize,
    const int maximumAdaptiveMeshDepth,
    const int haloCells,
    const int haloBufferCells,
    const int limiterBufferCells,
    const int regularisedFineGridLevels,
    const exahype::solvers::Solver::TimeStepping timeStepping,
    const int DMPObservables,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, Solver::Type::ADERDG, numberOfVariables,
             numberOfParameters, numberOfGlobalObservables, basisSize,
             maximumMeshSize, maximumAdaptiveMeshDepth,
             timeStepping, std::move(profiler)),
     _previousMinTimeStamp( std::numeric_limits<double>::infinity() ),
     _previousMinTimeStepSize( std::numeric_limits<double>::infinity() ),
     _minTimeStamp( std::numeric_limits<double>::infinity() ),
     _minTimeStepSize( std::numeric_limits<double>::infinity() ),
     _estimatedTimeStepSize( std::numeric_limits<double>::infinity() ),
     _admissibleTimeStepSize( std::numeric_limits<double>::infinity() ),
     _stabilityConditionWasViolated( false ),
     _minimumRefinementStatusToRequestMeshRefinementInVirtualCell(1+haloBufferCells),
     _refineOrKeepOnFineGrid(1+haloCells+haloBufferCells),
     _DMPObservables(DMPObservables),
     _minRefinementStatusForTroubledCell(_refineOrKeepOnFineGrid+3),
     _checkForNaNs(true),
      _meshUpdateEvent(MeshUpdateEvent::None) {
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

}

int exahype::solvers::ADERDGSolver::getUnknownsPerFace() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCellBoundary() const {
  return DIMENSIONS_TIMES_TWO * getUnknownsPerFace();
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getUnknownsPerCell(); // +1 for sources
}

int exahype::solvers::ADERDGSolver::getSpaceTimeUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getSpaceTimeUnknownsPerCell();  // +1 for sources
}

int exahype::solvers::ADERDGSolver::getDataPerFace() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getDataPerCellBoundary() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO;
}

int exahype::solvers::ADERDGSolver::getDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getDMPObservables() const {
  return _DMPObservables;
}

int exahype::solvers::ADERDGSolver::getMinRefinementStatusForTroubledCell() const {
  return _minRefinementStatusForTroubledCell;
}

void exahype::solvers::ADERDGSolver::resetMeshUpdateEvent() {
  _meshUpdateEvent = MeshUpdateEvent::None;
}

void exahype::solvers::ADERDGSolver::updateMeshUpdateEvent(MeshUpdateEvent meshUpdateEvent) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _meshUpdateEvent = mergeMeshUpdateEvents(_meshUpdateEvent,meshUpdateEvent);
  lock.free();
}

exahype::solvers::ADERDGSolver::MeshUpdateEvent
exahype::solvers::ADERDGSolver::getMeshUpdateEvent() const {
  return _meshUpdateEvent;
}

std::tuple<double,double> exahype::solvers::ADERDGSolver::getRiemannSolverTimeStepData(
    const CellDescription& cellDescription1,
    const CellDescription& cellDescription2) const {
  auto result  = std::make_tuple<double,double>(0.0,0.0);
  switch (_timeStepping) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      std::get<0>(result) = _minTimeStamp;
      std::get<1>(result) = _minTimeStepSize;
      break;
    default:
//      std::get<0>(result) = std::max(cellDescription1.getTimeStamp(),cellDescription2.getTimeStamp());
//      std::get<1>(result) = std::min(cellDescription1.getTimeStepSize(),cellDescription2.getTimeStepSize());
      logError("getRiemannSolverTimeStepData(...)","Unknown time stepping scheme.")
      std::abort();
  }
  return result;
}

std::tuple<double,double> exahype::solvers::ADERDGSolver::getPredictionTimeStepData(
    const CellDescription& cellDescription,const bool duringFusedTimeStep) const {
  auto result  = std::make_tuple<double,double>(0.0,0.0);
  switch (_timeStepping) {
    case TimeStepping::Global:
      if ( duringFusedTimeStep ) {
        std::get<0>(result) = _minTimeStamp+_minTimeStepSize;
        std::get<1>(result) = _estimatedTimeStepSize;
      } else {
        std::get<0>(result) = _minTimeStamp;
        std::get<1>(result) = _minTimeStepSize;
      }
      break;
    case TimeStepping::GlobalFixed:
      if ( duringFusedTimeStep ) {
        std::get<0>(result) = _minTimeStamp+_minTimeStepSize;
        std::get<1>(result) = _minTimeStepSize;
      } else {
        std::get<0>(result) = _minTimeStamp;
        std::get<1>(result) = _minTimeStepSize;
      }
      break;
      break;
    default:
      logError("getRiemannSolverTimeStepData(...)","Unknown time stepping scheme.")
      std::abort();
  }
  return result;
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
    CellDescription& p) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      p.setPreviousTimeStamp(_previousMinTimeStamp);
      p.setPreviousTimeStepSize(_previousMinTimeStepSize);
  
      p.setTimeStamp(_minTimeStamp);
      p.setTimeStepSize(_minTimeStepSize);
  }
}

void exahype::solvers::ADERDGSolver::kickOffTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _meshUpdateEvent               = MeshUpdateEvent::None;
    _admissibleTimeStepSize        = std::numeric_limits<double>::infinity();
    _stabilityConditionWasViolated = false;
    resetGlobalObservables(_nextGlobalObservables.data());
  }

  // call user code
  beginTimeStep(_minTimeStamp,isFirstTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::ADERDGSolver::wrapUpTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch,const bool isLastTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _previousMinTimeStepSize  = _minTimeStepSize;
    _previousMinTimeStamp     = _minTimeStamp;
  }
  _minTimeStamp += _minTimeStepSize;

  _stabilityConditionWasViolated = false;
  if (
      tarch::parallel::Node::getInstance().isGlobalMaster() &&
      getTimeStepping() != TimeStepping::GlobalFixed // fix the time step size in intermediate batch iterations
  ) {
    if ( FuseAllADERDGPhases && !isLinear() ) {
      if ( isLastTimeStepOfBatchOrNoBatch ) {
        if ( _estimatedTimeStepSize > _admissibleTimeStepSize ) { // rerun
          _minTimeStepSize       = FusedTimeSteppingRerunFactor * _admissibleTimeStepSize;
          _estimatedTimeStepSize = _minTimeStepSize;
          _stabilityConditionWasViolated = true;
        } else {
          _minTimeStepSize       = _estimatedTimeStepSize; // as we have computed the predictor with an estimate, we have to use the estimated time step size to perform the face integral
          _estimatedTimeStepSize = 0.5 * ( FusedTimeSteppingDiffusionFactor * _admissibleTimeStepSize + _estimatedTimeStepSize );
        }
      } else { // use fixed time step size in intermediate batch iterations
        _minTimeStepSize  = _estimatedTimeStepSize;
      }
    } else if ( !isLinear() ) { // non-fused, non-linear
      _minTimeStepSize = _admissibleTimeStepSize;
    } // else if linear do not change the time step size at all
  }

  if ( isLastTimeStepOfBatchOrNoBatch ) {
    std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),_globalObservables.begin());
    wrapUpGlobalObservables(_globalObservables.data());
  }

  // call user code
  endTimeStep(_minTimeStamp,isLastTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::ADERDGSolver::updateTimeStepSize() {
  if ( FuseAllADERDGPhases ) {
    _minTimeStepSize        = FusedTimeSteppingRerunFactor * _admissibleTimeStepSize;
    _estimatedTimeStepSize  = _minTimeStepSize;
  } else {
    _minTimeStepSize        = _admissibleTimeStepSize;
  }
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
  _stabilityConditionWasViolated = false;
}

void exahype::solvers::ADERDGSolver::updateGlobalObservables() {
  std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),
            _globalObservables.begin());
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep() {
  _minTimeStamp    = _previousMinTimeStamp;
  _minTimeStepSize = _previousMinTimeStepSize;

  _previousMinTimeStamp    = std::numeric_limits<double>::infinity();
  _previousMinTimeStepSize = std::numeric_limits<double>::infinity();

  // TODO(Lukas) Maybe also rollback global observables?
}

double exahype::solvers::ADERDGSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getEstimatedTimeStepSize() const {
  return _estimatedTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinTimeStepSize() const {
  return _previousMinTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinTimeStamp() const {
  return _previousMinTimeStamp;
}

double exahype::solvers::ADERDGSolver::getAdmissibleTimeStepSize() const {
  return _admissibleTimeStepSize;
}

void exahype::solvers::ADERDGSolver::updateAdmissibleTimeStepSize( double value ) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _admissibleTimeStepSize = std::min(_admissibleTimeStepSize,value);
  lock.free();
}

void exahype::solvers::ADERDGSolver::resetAdmissibleTimeStepSize() {
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
}

void exahype::solvers::ADERDGSolver::initSolver(
    const double                                timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const double                                boundingBoxSize,
    const double                                boundingBoxMeshSize,
    const std::vector<std::string>&             cmdlineargs,
    const exahype::parser::ParserView&          parserView
) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(
          std::min(boundingBoxMeshSize,_maximumMeshSize),boundingBoxSize);
  _coarsestMeshSize  = coarsestMeshInfo.first;
  _coarsestMeshLevel = coarsestMeshInfo.second;

  _previousMinTimeStamp = timeStamp;
  _minTimeStamp         = timeStamp;

  _previousMinTimeStepSize = 0.0;
  _minTimeStepSize         = 0.0;

  _estimatedTimeStepSize   = std::numeric_limits<double>::infinity();
  _admissibleTimeStepSize  = std::numeric_limits<double>::infinity();

  _meshUpdateEvent = MeshUpdateEvent::InitialRefinementRequested;

  // global observables
  _globalObservables.resize(_numberOfGlobalObservables);
  _nextGlobalObservables.resize(_numberOfGlobalObservables);
  resetGlobalObservables(_globalObservables.data()    );
  resetGlobalObservables(_nextGlobalObservables.data());

  init(cmdlineargs,parserView); // call user define initalisiation

  
  #ifdef Parallel
  _invalidExtrapolatedPredictor.resize(getBndFaceSize());
  _invalidFluctuations.resize(getBndFluxSize());
  std::fill_n(_invalidExtrapolatedPredictor.data(),_invalidExtrapolatedPredictor.size(),-1);
  std::fill_n(_invalidFluctuations.data(),_invalidFluctuations.size(),-1);

  _receivedExtrapolatedPredictor.resize(getBndFaceSize());
  _receivedFluctuations.resize(getBndFluxSize());

  _receivedUpdate.reserve(getUpdateSize());
  #endif
}

bool exahype::solvers::ADERDGSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  bool isPerformingPrediction = false;

  switch (section) {
    case exahype::State::AlgorithmSection::TimeStepping:
      isPerformingPrediction = true;
      break;
    case exahype::State::AlgorithmSection::PredictionRerunAllSend:
      isPerformingPrediction = !hasRequestedAnyMeshRefinement() &&
                               getStabilityConditionWasViolated();
      break;
    case exahype::State::AlgorithmSection::PredictionOrLocalRecomputationAllSend:
      isPerformingPrediction = hasRequestedAnyMeshRefinement();
      break;
    default:
      break;
  }

  return isPerformingPrediction;
}

bool exahype::solvers::ADERDGSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  return ( exahype::State::AlgorithmSection::MeshRefinement==section ) &&
         hasRequestedAnyMeshRefinement();
}

bool exahype::solvers::ADERDGSolver::getStabilityConditionWasViolated() const {
  return  _stabilityConditionWasViolated;
}

void exahype::solvers::ADERDGSolver::setStabilityConditionWasViolated(const bool state) {
  _stabilityConditionWasViolated = state;
}

bool exahype::solvers::ADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) {
  bool result = cellDescriptionsIndex>=0;
  assertion1(!result || Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  return result;
}

int exahype::solvers::ADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

///////////////////////////////////
// CELL-LOCAL MESH REFINEMENT
///////////////////////////////////


////////////////////////////////////////
// CELL-LOCAL
////////////////////////////////////////
void exahype::solvers::ADERDGSolver::validateCellDescriptionData(
  const CellDescription& cellDescription,
  const bool validateTimeStepData,
  const bool afterCompression,
  const bool beforePrediction,
  const std::string& methodTraceOfCaller) const {
  #ifdef Asserts
  if ( _checkForNaNs && validateTimeStepData ) {
    assertion2(std::isfinite(cellDescription.getTimeStamp()),
        cellDescription.toString(),toString());
    assertion2(cellDescription.getTimeStamp()>=0,
        cellDescription.toString(),toString());
  }
  // TODO(Dominic): Remove assertion1(cellDescription.getRefinementEvent()==CellDescription::None,cellDescription.toString());
  assertion1(getType()==exahype::solvers::Solver::Type::ADERDG,cellDescription.toString());

  if ( _checkForNaNs && afterCompression) {
    // TODO(Dominic)
  } else if ( _checkForNaNs ) {
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()),cellDescription.toString());

    double* luh  = static_cast<double*>(cellDescription.getSolution());
    double* lduh = static_cast<double*>(cellDescription.getUpdate());

    double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
    double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

    int dataPerCell             = getDataPerCell();
    int updateSize              = getUpdateSize();

    int dataPerCellBoundary     = getBndTotalSize();
    int unknownsPerCellBoundary = getBndFluxTotalSize();

    for (int i=0; i<dataPerCell; i++) {
      assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(luh[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }
    if ( !beforePrediction ) {
      for (int i=0; i<updateSize; i++) {
       assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lduh[i]),
           cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }

      for (int i=0; i<dataPerCellBoundary; i++) {
        assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lQhbnd[i]),
            cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }

      for (int i=0; i<unknownsPerCellBoundary; i++) {
        assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lFhbnd[i]),
            cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }
    }
  }
  #endif
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::ADERDGSolver::updateRefinementStatusAfterSolutionUpdate(CellDescription& cellDescription) {
  if ( OnlyInitialMeshRefinement ) {
    updateRefinementStatus(cellDescription);
    return MeshUpdateEvent::None;
  }

  cellDescription.setRefinementFlag(false);
  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    const double* solution = static_cast<double*>(cellDescription.getSolution());
    RefinementControl refinementControl = refinementCriterion(
                      solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
                      cellDescription.getSize(),
                      cellDescription.getTimeStamp(), // must be called after advancing in time
                      cellDescription.getLevel());
    if (
        (refinementControl==RefinementControl::Refine ||
        (cellDescription.getLevel()==getMaximumAdaptiveMeshLevel()    &&
        refinementControl==RefinementControl::Keep))
    ) {
      cellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
      cellDescription.setRefinementFlag(true);
    } else if (
        cellDescription.getRefinementStatus()<_refineOrKeepOnFineGrid &&
        cellDescription.getLevel()<getMaximumAdaptiveMeshLevel()     &&
        refinementControl==RefinementControl::Keep
    ) {
      cellDescription.setRefinementStatus(Keep);
    }
    else if (
        cellDescription.getRefinementStatus()<=Keep &&
        refinementControl==RefinementControl::Erase
    ) {
      cellDescription.setRefinementStatus(Erase);
    }

    // update refinement status after prescribing refinement values
    updateRefinementStatus(cellDescription);

    return
        (cellDescription.getLevel() < getMaximumAdaptiveMeshLevel() &&
         refinementControl==RefinementControl::Refine ) ?
            MeshUpdateEvent::RefinementRequested : MeshUpdateEvent::None;
  } else if ( cellDescription.getType()==CellDescription::Type::Virtual ) {
    // bottom up refinement criterion TODO(Dominic): Add to docu
    // We allow the halo region to diffuse into the virtual subcells
    // up to some point.
    updateRefinementStatus(cellDescription);
    if (
        cellDescription.getLevel()==getMaximumAdaptiveMeshLevel() &&
        cellDescription.getRefinementStatus() >= _minimumRefinementStatusToRequestMeshRefinementInVirtualCell
    ) {
      return MeshUpdateEvent::RefinementRequested;
    } else {
      return MeshUpdateEvent::None;
    }
  } else {
    cellDescription.setRefinementStatus(Erase); // Cannot override the refinement / limiter status in other cells
    return MeshUpdateEvent::None;
  }
}

void exahype::solvers::ADERDGSolver::fusedTimeStepBody(
    CellDescription&                                   cellDescription,
    CellInfo&                                          cellInfo,
    const double                                       predictionTimeStamp,
    const double                                       predictionTimeStepSize,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         isSkeletonCell,
    const bool                                         mustBeDoneImmediately) {
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_begin(fusedTimeStepBodyHandleSkeleton);
  } else {
    VT_begin(fusedTimeStepBodyHandle);
  }
  #endif

  correction(cellDescription,boundaryMarkers,isFirstTimeStepOfBatch,isFirstTimeStepOfBatch/*addSurfaceIntegralContributionToUpdate*/);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(cellDescription,isFirstTimeStepOfBatch);
  cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(cellDescription);

  reduce(cellDescription,result);

  if (
      SpawnPredictionAsBackgroundJob &&
      !mustBeDoneImmediately &&
      isLastTimeStepOfBatch // only spawned in last iteration if a FusedTimeStepJob was spawned before
  ) {
    const int element = cellInfo.indexOfADERDGCellDescription(cellDescription.getSolverNumber());
    peano::datatraversal::TaskSet( new PredictionJob(
        *this, cellDescription, cellInfo._cellDescriptionsIndex,element,
        predictionTimeStamp, predictionTimeStepSize,
        false/*is uncompressed*/, isSkeletonCell, isLastTimeStepOfBatch/*addVolumeIntegralResultToUpdate*/ ) );
  } else {
    predictionAndVolumeIntegralBody(
        cellDescription,
        predictionTimeStamp, predictionTimeStepSize,
        false, isSkeletonCell, isLastTimeStepOfBatch/*addVolumeIntegralResultToUpdate*/);
  }

  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_end(fusedTimeStepBodyHandleSkeleton);
  } else {
    VT_end(fusedTimeStepBodyHandle);
  }
  #endif
}

void exahype::solvers::ADERDGSolver::fusedTimeStepOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      const bool isAMRSkeletonCell     = belongsToAMRSkeleton(cellDescription);
      const bool isAtRemoteBoundary    = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
      const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
      const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

      if (
          (SpawnUpdateAsBackgroundJob || (SpawnPredictionAsBackgroundJob && !isLastTimeStepOfBatch)) &&
          !mustBeDoneImmediately
      ) {
        const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,true);
        peano::datatraversal::TaskSet( new FusedTimeStepJob(
            *this, cellDescription, cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch, isLastTimeStepOfBatch,
            boundaryMarkers, isSkeletonCell) );
      } else {
        const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,true);
        fusedTimeStepBody(
            cellDescription,cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch,isLastTimeStepOfBatch,
            boundaryMarkers,isSkeletonCell,mustBeDoneImmediately );
      }
    } else if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication
    ) {
      restrictToTopMostParent(cellDescription,isFirstTimeStepOfBatch/*addToCoarseGridUpdate*/);
      updateMeshUpdateEvent( updateRefinementStatusAfterSolutionUpdate(cellDescription) );
      cellDescription.setHasCompletedLastStep(true);
    } else {
      cellDescription.setHasCompletedLastStep(true);
    }
  }
}

void exahype::solvers::ADERDGSolver::reduce(
    const CellDescription& cellDescription,
    const UpdateResult&    result) {
  updateMeshUpdateEvent(result._meshUpdateEvent);
  updateAdmissibleTimeStepSize(result._timeStepSize);

  if ( _numberOfGlobalObservables > 0 ) {
    assert(cellDescription.getType()==CellDescription::Type::Leaf);
    const double* const luh = static_cast<double*>(cellDescription.getSolution());
    const auto cellCentre   = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
    const auto& cellSize    = cellDescription.getSize();
    const auto t             = cellDescription.getTimeStamp();
    const auto dt            = cellDescription.getTimeStepSize();
    updateGlobalObservables(_nextGlobalObservables.data(),luh,cellCentre,cellSize,t,dt);
  }
}

void exahype::solvers::ADERDGSolver::updateBody(
    CellDescription&                                   cellDescription,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  #ifdef USE_ITAC
  VT_begin(updateBodyHandle);
  #endif

  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  uncompress(cellDescription);

  correction(cellDescription,boundaryMarkers,true,false/*effect: add face integral result directly to solution*/);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(cellDescription,true);
  cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(cellDescription);

  reduce(cellDescription,result);

  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
  compress(cellDescription,isAtRemoteBoundary);

  cellDescription.setHasCompletedLastStep(true); // required as prediction checks the flag too. Field should be renamed "setHasCompletedLastOperation(...)".

  #ifdef USE_ITAC
  VT_end(updateBodyHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::updateOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
    if ( cellDescription.getType()==CellDescription::Type::Leaf && SpawnUpdateAsBackgroundJob ) {
      peano::datatraversal::TaskSet ( new UpdateJob(*this,cellDescription,cellInfo,boundaryMarkers) );
    }
    else if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      updateBody(cellDescription,cellInfo,boundaryMarkers);
    }
    else if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication
    ) {
      restrictToTopMostParent(cellDescription,false/*effect: add face integral result directly to solution*/);
      updateMeshUpdateEvent( updateRefinementStatusAfterSolutionUpdate(cellDescription) );
      cellDescription.setHasCompletedLastStep(true);
    }
    else {
      cellDescription.setHasCompletedLastStep(true);
    }
  }
}

void exahype::solvers::ADERDGSolver::compress(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      const bool isSkeletonCell = belongsToAMRSkeleton(cellDescription);
      compress(cellDescription,isSkeletonCell);
    }
  }
}

int exahype::solvers::ADERDGSolver::predictionAndVolumeIntegralBody(
    CellDescription& cellDescription,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isSkeletonCell,
    const bool   addVolumeIntegralResultToUpdate) {
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_begin(predictorBodyHandleSkeleton);
  } else {
    VT_begin(predictorBodyHandle);
  }
  #endif
  if (uncompressBefore) { uncompress(cellDescription); }

  validateCellDescriptionData(cellDescription,true,false,true,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [pre]");

  double* luh    = static_cast<double*>(cellDescription.getSolution());
  double* lduh   = static_cast<double*>(cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("performPredictionAndVolumeIntegralBody(...)","#predictions="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  const int numberOfPicardIterations = fusedSpaceTimePredictorVolumeIntegral(
      lduh,lQhbnd,lGradQhbnd,lFhbnd,
      luh,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      predictorTimeStamp,
      predictorTimeStepSize,
      addVolumeIntegralResultToUpdate); // TODO(Dominic): fix 'false' case

  compress(cellDescription,isSkeletonCell);

  validateCellDescriptionData(cellDescription,true,true,false,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [post]");

  cellDescription.setHasCompletedLastStep(true);

  
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_end(predictorBodyHandleSkeleton);
  } else {
    VT_end(predictorBodyHandle);
  }
  #endif

  return numberOfPicardIterations;
}

void exahype::solvers::ADERDGSolver::predictionAndVolumeIntegral(
    const int    solverNumber,
    CellInfo&    cellInfo,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isAtRemoteBoundary,
    const bool   addVolumeIntegralResultToUpdate) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    cellDescription.setHasCompletedLastStep(false);

    const bool isAMRSkeletonCell     = belongsToAMRSkeleton(cellDescription);
    const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
    const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

    if ( SpawnPredictionAsBackgroundJob && !mustBeDoneImmediately ) {
      peano::datatraversal::TaskSet( new PredictionJob(
              *this, cellDescription, cellInfo._cellDescriptionsIndex, element,
              predictorTimeStamp,predictorTimeStepSize,
              uncompressBefore,isSkeletonCell,addVolumeIntegralResultToUpdate) );
    }
    else {
      predictionAndVolumeIntegralBody(
          cellDescription,
          predictorTimeStamp,predictorTimeStepSize,
          uncompressBefore,isSkeletonCell,addVolumeIntegralResultToUpdate);
    }
  }
}

void exahype::solvers::ADERDGSolver::predictionAndVolumeIntegral(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);

    const bool isAMRSkeletonCell = belongsToAMRSkeleton(cellDescription);
    const bool isSkeletonCell    = isAMRSkeletonCell || isAtRemoteBoundary;
    waitUntilCompletedLastStep(cellDescription,isSkeletonCell,false);
    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,false); // this is either the fused scheme or a predictor recomputation

      const bool rollbacksPossible = !OnlyInitialMeshRefinement;
      if ( !FuseAllADERDGPhases && rollbacksPossible ) { // backup previous solution here as prediction already adds a contribution to solution if not all alg. phases are fused.
        std::copy_n(
            static_cast<double*>(cellDescription.getSolution()),getDataPerCell(),
            static_cast<double*>(cellDescription.getPreviousSolution()));
      }

      predictionAndVolumeIntegral(solverNumber,cellInfo,
          std::get<0>(predictionTimeStepData),
          std::get<1>(predictionTimeStepData),
          true,isAtRemoteBoundary,
          FuseAllADERDGPhases/*addVolumeIntegralResultToUpdate*/);
    }
  }
}

double exahype::solvers::ADERDGSolver::computeTimeStepSize(CellDescription& cellDescription) {
  if( cellDescription.getType()==CellDescription::Type::Leaf ) {
    const double* luh = static_cast<double*>(cellDescription.getSolution());

    validateCellDescriptionData(cellDescription,false,false,true,"computeTimeStepSizes(...)");
    double admissibleTimeStepSize = stableTimeStepSize(luh,cellDescription.getSize());
    assertion2(!_checkForNaNs || admissibleTimeStepSize>0,admissibleTimeStepSize,cellDescription.toString());

    assertion3(!_checkForNaNs || admissibleTimeStepSize<std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),admissibleTimeStepSize,cellDescription.toString());
    assertion2(!_checkForNaNs || std::isfinite(admissibleTimeStepSize),admissibleTimeStepSize,cellDescription.toString());

    return admissibleTimeStepSize;
  } else {
    return std::numeric_limits<double>::infinity();
  }
}


double exahype::solvers::ADERDGSolver::startNewTimeStep(
    CellDescription& cellDescription,
    const bool isFirstTimeStepOfBatch) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  // n-1
  if ( isFirstTimeStepOfBatch ) {
    cellDescription.setPreviousTimeStamp(cellDescription.getTimeStamp());
    cellDescription.setPreviousTimeStepSize(cellDescription.getTimeStepSize());
  }
  // n
  cellDescription.setTimeStamp(cellDescription.getTimeStamp()+cellDescription.getTimeStepSize());
  if ( getTimeStepping() != TimeStepping::GlobalFixed ) {
    double admissibleTimeStepSize = computeTimeStepSize(cellDescription);
    cellDescription.setTimeStepSize(admissibleTimeStepSize);
    return admissibleTimeStepSize;
  } else {
    return cellDescription.getTimeStepSize();
  }
}

void exahype::solvers::ADERDGSolver::updateTimeStepSize(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    const double admissibleTimeStepSize = computeTimeStepSize(cellDescription);
    cellDescription.setTimeStepSize( admissibleTimeStepSize );
    cellDescription.setHasCompletedLastStep(true);
    updateAdmissibleTimeStepSize(admissibleTimeStepSize);
  }
}

void exahype::solvers::ADERDGSolver::updateGlobalObservables(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if ( _numberOfGlobalObservables > 0 && cellDescription.getType() == CellDescription::Type::Leaf ) {
      const double* const luh  = static_cast<double*>(cellDescription.getSolution());
      const auto cellCentre    = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
      const auto& cellSize     = cellDescription.getSize();
      const auto t             = cellDescription.getTimeStamp();
      const auto dt            = cellDescription.getTimeStepSize();
      updateGlobalObservables(_nextGlobalObservables.data(),luh,cellCentre,cellSize,t,dt);
    }
  }
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep(CellDescription& cellDescription) const {
  cellDescription.setTimeStamp   (cellDescription.getPreviousTimeStamp());
  cellDescription.setTimeStepSize(std::numeric_limits<double>::infinity());

  cellDescription.setPreviousTimeStamp(std::numeric_limits<double>::infinity());
  cellDescription.setPreviousTimeStepSize(std::numeric_limits<double>::infinity());
}

void exahype::solvers::ADERDGSolver::adjustSolution(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf ||
             cellDescription.getType()==CellDescription::Type::LeafProlongates ||
             cellDescription.getType()==CellDescription::Type::ParentCoarsens
             ,cellDescription.toString());

  double* solution = static_cast<double*>(cellDescription.getSolution());
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp(),
      cellDescription.getTimeStepSize());

  double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
  adjustSolution(
      previousSolution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getPreviousTimeStamp(),
      cellDescription.getPreviousTimeStepSize());

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getDataPerCell(); i++) {
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"adjustSolution(...)",i);
    }
  }
  #endif
}


void exahype::solvers::ADERDGSolver::printADERDGSolution2D(const CellDescription& cellDescription)  const {
  #if DIMENSIONS==2
  double* solution = static_cast<double*>(cellDescription.getSolution());
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  assertion1(solution!=nullptr,cellDescription.toString());

  std::cout <<  "solution:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int unknown=0; unknown < numberOfData; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis)*numberOfData+unknown;
      std::cout << std::setprecision(3) << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::ADERDGSolver::printADERDGExtrapolatedPredictor2D(const CellDescription& cellDescription) const {
  #if DIMENSIONS==2
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  assertion1(extrapolatedPredictor!=nullptr,cellDescription.toString());

  std::cout <<  "extrapolated predictor:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int f=0; f<DIMENSIONS_TIMES_TWO; f++) {
    std::cout <<  "face=" << f << std::endl;
    for (int unknown=0; unknown < numberOfData; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for (int i=0; i<_nodesPerCoordinateAxis; i++) {
        int iScalar = f*numberOfData*_nodesPerCoordinateAxis + numberOfData*i + unknown;
        std::cout << std::setprecision(3) << extrapolatedPredictor[iScalar] << ",";
        if (i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::ADERDGSolver::printADERDGFluctuations2D(const CellDescription& cellDescription) const {
  #if DIMENSIONS==2
  double* fluctuation = static_cast<double*>(cellDescription.getFluctuation());
  assertion1(fluctuation!=nullptr,cellDescription.toString());

  std::cout << "fluctuations:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int f=0; f<DIMENSIONS_TIMES_TWO; f++) {
    std::cout <<  "face=" << f << std::endl;
    for (int unknown=0; unknown < numberOfData; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for (int i=0; i<_nodesPerCoordinateAxis; i++) {
        int iScalar = f*numberOfData*_nodesPerCoordinateAxis + numberOfData*i + unknown;
        std::cout << std::setprecision(3) << fluctuation[iScalar] << ",";
        if (i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}



void exahype::solvers::ADERDGSolver::surfaceIntegral(
    CellDescription&                                   cellDescription,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         addToUpdate) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  double* output= static_cast<double*>(cellDescription.getSolution());
  if ( addToUpdate ) { output = static_cast<double*>(cellDescription.getUpdate()); }

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(tarch::la::equals(cellDescription.getTimeStepSize(),0.0)  || std::isfinite(output[i]),cellDescription.toString(),"surfaceIntegral",i);
    }
  }
  #endif

  // the actual operation
  const int dofsPerFace = getBndFluxSize();
  for (int direction=0; direction<DIMENSIONS; direction++) {
    for (int orientation=0; orientation<2; orientation++) {
      const int faceIndex=2*direction+orientation;
      // impose boundary conditions
      if ( boundaryMarkers[faceIndex]==mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex ) {
        mergeWithBoundaryData(cellDescription,faceIndex,direction,orientation);
      }
      // perform face integral
      if ( cellDescription.getFacewiseAugmentationStatus(faceIndex)<MaximumAugmentationStatus ) { // ignore Ancestors
        double* const lFhbnd = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * faceIndex;
        faceIntegral(output,lFhbnd,direction,orientation,0/*implicit conversion*/,0,cellDescription.getSize(),
                     cellDescription.getTimeStepSize(),addToUpdate);
      }
    }
  }

  // check that all neighbour merges have been performed; reset the flags
  if ( !tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) && !ProfileUpdate ) {
    logError("surfaceIntegral(...)","Riemann solve was not performed on all faces of cell= "<<cellDescription.toString());
    std::terminate();
  }
  assertion1( tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) || ProfileUpdate,cellDescription.toString());
  cellDescription.setNeighbourMergePerformed(static_cast<signed char>(false));

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(std::isfinite(output[i]),cellDescription.toString(),"updateSolution(...)",i);
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::addUpdateToSolution(CellDescription& cellDescription,const bool backupPreviousSolution) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  double* update = static_cast<double*>(cellDescription.getUpdate());

  if ( backupPreviousSolution ) {
    // perform the update
    swapSolutionAndPreviousSolution(cellDescription); // solution is overwritten with the current solution plus the update,
    // while current solution is remembered as the previous solution.
    double* solution         = static_cast<double*>(cellDescription.getSolution());
    double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
    addUpdateToSolution(solution,previousSolution,update,cellDescription.getTimeStepSize());
  } else {
    double* solution = static_cast<double*>(cellDescription.getSolution());
    addUpdateToSolution(solution,solution,update,cellDescription.getTimeStepSize());
  }
}

void exahype::solvers::ADERDGSolver::adjustSolutionAfterUpdate(CellDescription& cellDescription) {
  double* solution = static_cast<double*>(cellDescription.getSolution());
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),
      cellDescription.getTimeStepSize());

  // only for profiling
  if ( Solver::ProfileUpdate ) { swapSolutionAndPreviousSolution(cellDescription); }

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"updateSolution(...)",i);
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::correction(
    CellDescription&                                   cellDescription,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         backupPreviousSolution,
    const bool                                         addSurfaceIntegralResultToUpdate) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString())
  assertion1(std::isfinite(cellDescription.getTimeStamp()   ),cellDescription.toString());
  assertion1(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString());
  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("mergeNeighboursData(...)","#updateSolution="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  surfaceIntegral(cellDescription,boundaryMarkers,addSurfaceIntegralResultToUpdate);
  if ( addSurfaceIntegralResultToUpdate ) {
    addUpdateToSolution(cellDescription,backupPreviousSolution);
  }
  adjustSolutionAfterUpdate(cellDescription);

  // update grid flags
  updateCommunicationStatus(cellDescription);
  updateAugmentationStatus(cellDescription);
}

void exahype::solvers::ADERDGSolver::swapSolutionAndPreviousSolution(CellDescription& cellDescription) const {
  assertion(cellDescription.getType()==CellDescription::Type::Leaf);

  // Simply swap the heap indices
  const int previousSolutionIndex = cellDescription.getPreviousSolutionIndex();
  void* previousSolution          = cellDescription.getPreviousSolution(); // pointer
  cellDescription.setPreviousSolutionIndex(cellDescription.getSolutionIndex());
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolutionIndex(previousSolutionIndex);
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::ADERDGSolver::prolongateObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& parentCellDescription) const {
  const int numberOfObservables = getDMPObservables();
  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; ++faceIndex) {
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) { // TODO(Dominic): If the grid changes dynamically during the time steps,
      // fine
      double* minFine = static_cast<double*>(cellDescription.getSolutionMin()) + numberOfObservables * faceIndex;
      double* maxFine = static_cast<double*>(cellDescription.getSolutionMax()) + numberOfObservables * faceIndex;
      // coarse
      const double* minCoarse = static_cast<double*>(parentCellDescription.getSolutionMin()) +  numberOfObservables * faceIndex;
      const double* maxCoarse = static_cast<double*>(parentCellDescription.getSolutionMax()) +  numberOfObservables * faceIndex;

      std::copy_n( minCoarse,numberOfObservables, minFine );
      std::copy_n( maxCoarse,numberOfObservables, maxFine );
    }
  }
}

void exahype::solvers::ADERDGSolver::prolongateFaceDataToVirtualCell(
    CellDescription& cellDescription,
    const CellDescription& parentCellDescription,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  #ifdef USE_ITAC
  VT_begin(prolongateFaceDataToVirtualCellHandle);
  #endif

  assertion(parentCellDescription.getSolverNumber() == cellDescription.getSolverNumber());
  assertion(parentCellDescription.getType() == CellDescription::Type::Leaf ||
            parentCellDescription.getType() == CellDescription::Type::Virtual);

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = parentCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  double* update = static_cast<double*>(cellDescription.getUpdate());
  std::fill_n(update,getUpdateSize(),0.0);

  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; ++faceIndex) {
    const int direction = faceIndex/2;
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) { // TODO(Dominic): If the grid changes dynamically during the time steps,
      #ifdef Parallel
      assertion4( exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelFine-levelCoarse),
            cellDescription.toString(),parentCellDescription.toString(),tarch::parallel::Node::getInstance().getRank(),
            tarch::parallel::NodePool::getInstance().getMasterRank() ); // necessary but not sufficient
      #else
      assertion2( exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelFine-levelCoarse),
                  cellDescription.toString(),parentCellDescription.toString() ); // necessary but not sufficient
      #endif

      logDebug("prolongateFaceDataToDescendant(...)","cell=" << cellDescription.getOffset() <<
               ",level=" << cellDescription.getLevel() <<
               ",face=" << faceIndex <<
               ",subcellIndex" << subcellIndex.toString() <<
               " from " <<
               " parentCell="<<parentCellDescription.getOffset()<<
               " level="<<parentCellDescription.getLevel());

      // extrapolated predictor and flux interpolation
      // extrapolated predictor
      assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
      assertion1(DataHeap::getInstance().isValidIndex(parentCellDescription.getExtrapolatedPredictorIndex()),parentCellDescription.toString());

      const int dataPerFace = getBndFaceSize();
      const int dofsPerFace = getBndFluxSize();

      // fine
      double* lQhbndFine = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * faceIndex; // TODO(Dominic ): Pointer should be obtained only once
      double* lFhbndFine = static_cast<double*>(cellDescription.getFluctuation())           + dofsPerFace  * faceIndex ;
      // coarse
      const double* lQhbndCoarse = static_cast<double*>(parentCellDescription.getExtrapolatedPredictor()) + dataPerFace * faceIndex;
      const double* lFhbndCoarse = static_cast<double*>(parentCellDescription.getFluctuation()          ) + dofsPerFace  * faceIndex;

      faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,lFhbndCoarse, levelCoarse, levelFine,
                               exahype::amr::getSubfaceIndex(subcellIndex,direction));
    }
  }

  if ( getDMPObservables()>0 ) {
    prolongateObservablesMinAndMax(cellDescription,parentCellDescription);
  }

  cellDescription.setHasCompletedLastStep(true);

  #ifdef USE_ITAC
  VT_end(prolongateFaceDataToVirtualCellHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::prolongateFaceData(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);

  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication &&
        isValidCellDescriptionIndex(cellDescription.getParentIndex()) // might be at master-worker boundary
    ) {
        Solver::SubcellPosition subcellPosition = amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(cellDescription);
        CellDescription& parentCellDescription = getCellDescription(
            subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
        assertion1(parentCellDescription.getType()==CellDescription::Type::Leaf,parentCellDescription.toString());

        waitUntilCompletedLastStep<CellDescription>(parentCellDescription,true,false); // TODO(Dominic): We wait for skeleton jobs here. It might make sense to receiveDanglingMessages here too
        if (
            !SpawnProlongationAsBackgroundJob ||
            isAtRemoteBoundary
        ) {
          prolongateFaceDataToVirtualCell(cellDescription,parentCellDescription,subcellPosition.subcellIndex);
        } else {
          cellDescription.setHasCompletedLastStep(false); // done here in order to skip lookup of cell description in job constructor
          peano::datatraversal::TaskSet spawn( new ProlongationJob( *this, 
              cellDescription, parentCellDescription, subcellPosition.subcellIndex) );
        }
      }
      assertion2(
          cellDescription.getType()!=CellDescription::Type::Virtual ||
          isValidCellDescriptionIndex(cellDescription.getParentIndex()),
          cellDescription.toString(),
          tarch::parallel::Node::getInstance().getRank());
  }
}

void exahype::solvers::ADERDGSolver::restrictObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& parentCellDescription) const {
  const int numberOfObservables = getDMPObservables();
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) {
      // fine
      const double* minFine = static_cast<double*>(cellDescription.getSolutionMin()) + numberOfObservables* faceIndex;
      const double* maxFine = static_cast<double*>(cellDescription.getSolutionMax()) + numberOfObservables* faceIndex;
      // coarse
      double* minCoarse = static_cast<double*>(parentCellDescription.getSolutionMin()) + numberOfObservables * faceIndex;
      double* maxCoarse = static_cast<double*>(parentCellDescription.getSolutionMax()) + numberOfObservables * faceIndex;

      tarch::multicore::Lock lock(RestrictionSemaphore);
      for (int i=0; i<numberOfObservables; i++) {
        *(minCoarse+i) = std::min( *(minFine+i), *(minCoarse+i) );
        *(maxCoarse+i) = std::max( *(maxFine+i), *(maxCoarse+i) );
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::restrictToTopMostParent(
    CellDescription& cellDescription,
    const bool       addToCoarseGridUpdate) {
  #ifdef USE_ITAC
  VT_begin(restrictToTopMostParentHandle);
  #endif

  // validate and obtain parent
  assertion1( cellDescription.getType()==CellDescription::Type::Virtual &&
              cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication, cellDescription.toString() );
  assertion1( tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber()) != NotFound, cellDescription.toString());
  exahype::solvers::Solver::SubcellPosition subcellPosition =
      exahype::amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(cellDescription);
  assertion1(subcellPosition.parentElement!=exahype::solvers::Solver::NotFound,cellDescription.toString());

  CellDescription& parentCellDescription =
      getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);

  assertion(parentCellDescription.getSolverNumber()==cellDescription.getSolverNumber());
  assertion1(cellDescription.getType()==CellDescription::Type::Virtual,cellDescription.toString());
  assertion1(parentCellDescription.getType()==CellDescription::Type::Leaf,parentCellDescription.toString());


  // Perform the face integrals on fine grid cell
  double* updateFine = static_cast<double*>(cellDescription.getUpdate()); // TODO(Dominic): Can be temporary
  std::fill_n(updateFine,getUpdateSize(),0.0);

  const int levelDelta = cellDescription.getLevel() - parentCellDescription.getLevel();
  const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
      exahype::amr::computeSubcellIndex(
          cellDescription.getOffset(),cellDescription.getSize(),
          parentCellDescription.getOffset()); // TODO(Dominic): Maybe, I get can get rid of some variables again

  // gather contributions
  double dt = parentCellDescription.getTimeStepSize();
  switch (getTimeStepping()) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      // do nothing
      break;
    default:
      logError("restrictToTopMostParent(...)","Time stepping scheme not supported.");
      break;
  }

  const int dofsPerFace = getBndFluxSize();
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    const int direction   = faceIndex / 2;
    const int orientation = faceIndex % 2;
    const tarch::la::Vector<DIMENSIONS-1,int> subfaceIndex =
      exahype::amr::getSubfaceIndex(subcellIndex,direction);
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) {
      assertion1(exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelDelta),cellDescription.toString());
      assertion1(SpawnProlongationAsBackgroundJob || cellDescription.getNeighbourMergePerformed(faceIndex),cellDescription.toString());// necessary but not sufficient

      double* const lFhbnd = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * faceIndex;
      faceIntegral(updateFine,lFhbnd,direction,orientation,subfaceIndex,levelDelta,cellDescription.getSize(),dt,addToCoarseGridUpdate);

      logDebug("restrictToTopMostParent(...)","cell=" << cellDescription.getOffset() <<
             ",level=" << cellDescription.getLevel() <<
             ",face=" << faceIndex <<
             ",subcellIndex" << subcellIndex.toString() <<
             " to " <<
             " parentCell="<<parentCellDescription.getOffset()<<
             " level="<<parentCellDescription.getLevel());
    }
  }
  cellDescription.setNeighbourMergePerformed(static_cast<signed char>(false));

  // Add child contributions to parent
  if ( addToCoarseGridUpdate ) {
    double* updateCoarse = static_cast<double*>(parentCellDescription.getUpdate());
    tarch::multicore::Lock lock(RestrictionSemaphore);
    for (int i = 0; i < getUpdateSize(); ++i) { // TODO(Dominic): There must be a correction somewhere here or to the face integral if LTS is performed
      updateCoarse[i] += updateFine[i];
    }
    lock.free();
  } else {
    double* solutionCoarse = static_cast<double*>(parentCellDescription.getSolution());
    tarch::multicore::Lock lock(RestrictionSemaphore);
    addUpdateToSolution(solutionCoarse,solutionCoarse,updateFine,dt);
    lock.free();
  }

  if ( getDMPObservables()>0 ) {
    restrictObservablesMinAndMax(cellDescription,parentCellDescription);
  }

  #ifdef USE_ITAC
  VT_end(restrictToTopMostParentHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::disableCheckForNaNs() {
  _checkForNaNs = false;
}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::ADERDGSolver::rollbackSolutionGlobally(const int solverNumber,CellInfo& cellInfo) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    // 1. Rollback time step data
    rollbackToPreviousTimeStep(cellDescription);
    // 2. Rollback solution to previous one
    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      swapSolutionAndPreviousSolution(cellDescription);
    }

    // 3. Reset the previous refinement status
    cellDescription.setRefinementStatus(cellDescription.getPreviousRefinementStatus());
  }
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo, // corresponds to dest
    const int                                    neighbourAugmentationStatus,
    const int                                    neighbourCommunicationStatus,
    const int                                    neighbourRefinementStatus,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    const tarch::la::Vector<DIMENSIONS, int>&    posNeighbour,
    const tarch::la::Vector<DIMENSIONS, double>& barycentreFromVertex) {
  assertion(tarch::la::countEqualEntries(posNeighbour,pos)==DIMENSIONS-1); // only consider faces
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=Solver::NotFound ) {
    Solver::BoundaryFaceInfo face(pos,posNeighbour); //
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    const tarch::la::Vector<DIMENSIONS,double> barycentreFrom1 =
        exahype::Cell::computeFaceBarycentre(cellDescription.getOffset(),cellDescription.getSize(),face._direction,face._orientation);
    if ( Vertex::equalUpToRelativeTolerance(barycentreFrom1,barycentreFromVertex) ) {
      logDebug("mergeWithNeighbourMetadata(...)", "received neighbour metadata="<<neighbourAugmentationStatus<<","<<neighbourCommunicationStatus<<","<<neighbourRefinementStatus);

      mergeWithAugmentationStatus (cellDescription,face._faceIndex,neighbourAugmentationStatus );
      mergeWithCommunicationStatus(cellDescription,face._faceIndex,neighbourCommunicationStatus);
      mergeWithRefinementStatus   (cellDescription,face._faceIndex,neighbourRefinementStatus   );

      cellDescription.setNeighbourMergePerformed(face._faceIndex,true);
    }
  }
}

void exahype::solvers::ADERDGSolver::mergeNeighboursMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo1,
    Solver::CellInfo&                            cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>&    pos1,
    const tarch::la::Vector<DIMENSIONS, int>&    pos2,
    const tarch::la::Vector<DIMENSIONS, double>& barycentreFromVertex,
    const bool                                   checkThoroughly) const {
  const int element1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    Solver::InterfaceInfo face(pos1,pos2);
    CellDescription& cellDescription1 = cellInfo1._ADERDGCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._ADERDGCellDescriptions[element2];

    bool mergeMetadata = true;
    if ( checkThoroughly ) {
      const tarch::la::Vector<DIMENSIONS,double> barycentreFrom1 =
          exahype::Cell::computeFaceBarycentre(cellDescription1.getOffset(),cellDescription1.getSize(),face._direction,face._orientation1);
      const tarch::la::Vector<DIMENSIONS,double> barycentreFrom2 =
          exahype::Cell::computeFaceBarycentre(cellDescription2.getOffset(),cellDescription2.getSize(),face._direction,face._orientation2);
      mergeMetadata &= Vertex::equalUpToRelativeTolerance(barycentreFrom1,barycentreFromVertex) &&
                       Vertex::equalUpToRelativeTolerance(barycentreFrom2,barycentreFromVertex);
    }
    if ( mergeMetadata ) {
      mergeWithCommunicationStatus(cellDescription1,face._faceIndex1,cellDescription2.getCommunicationStatus());
      mergeWithAugmentationStatus( cellDescription1,face._faceIndex1,cellDescription2.getAugmentationStatus());
      mergeWithRefinementStatus(   cellDescription1,face._faceIndex1,cellDescription2.getRefinementStatus());

      mergeWithCommunicationStatus(cellDescription2,face._faceIndex2,cellDescription1.getCommunicationStatus());
      mergeWithAugmentationStatus( cellDescription2,face._faceIndex2,cellDescription1.getAugmentationStatus());
      mergeWithRefinementStatus(   cellDescription2,face._faceIndex2,cellDescription1.getRefinementStatus());
    } else {
      mergeWithCommunicationStatus(cellDescription1,face._faceIndex1,EmptyStatus);
      mergeWithAugmentationStatus( cellDescription1,face._faceIndex1,EmptyStatus);
      mergeWithRefinementStatus(   cellDescription1,face._faceIndex1,EmptyStatus);

      mergeWithCommunicationStatus(cellDescription2,face._faceIndex2,EmptyStatus);
      mergeWithAugmentationStatus( cellDescription2,face._faceIndex2,EmptyStatus);
      mergeWithRefinementStatus(   cellDescription2,face._faceIndex2,EmptyStatus);
    }

    cellDescription1.setNeighbourMergePerformed(face._faceIndex1,true); // here we only set, doesn't matter if operation is done twice.
    cellDescription2.setNeighbourMergePerformed(face._faceIndex2,true);
  }
}

// merge compute data
void exahype::solvers::ADERDGSolver::mergeNeighboursData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo1,
    Solver::CellInfo&                         cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  #ifdef USE_ITAC
  VT_begin(mergeNeighboursHandle);
  #endif

  const int element1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    Solver::InterfaceInfo face(pos1,pos2);
    CellDescription& cellDescription1 = cellInfo1._ADERDGCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._ADERDGCellDescriptions[element2];

    if ( ADERDGSolver::communicateWithNeighbour(cellDescription1,face._faceIndex1) ) {
      assertion1( ADERDGSolver::communicateWithNeighbour(cellDescription2,face._faceIndex2),cellDescription2.toString() );

      prefetchFaceData(cellDescription1,face._faceIndex1);
      prefetchFaceData(cellDescription2,face._faceIndex2);

      #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
      static int counter = 0;
      static double timeStamp = 0;
      if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
        logInfo("mergeNeighboursData(...)","#riemanns="<<counter);
        timeStamp = _minTimeStamp;
        counter=0;
      }
      counter++;
      #endif

      waitUntilCompletedLastStep<CellDescription>(cellDescription1,false,false);  // must be done before any other operation on the patches
      waitUntilCompletedLastStep<CellDescription>(cellDescription2,false,false);

      if ( CompressionAccuracy > 0.0 ) {
        peano::datatraversal::TaskSet uncompression(
            [&] () -> bool {
          uncompress(cellDescription1);
          return false;
        },
        [&] () -> bool {
          uncompress(cellDescription1);
          return false;
        },
        peano::datatraversal::TaskSet::TaskType::Background,
        peano::datatraversal::TaskSet::TaskType::Background,
        true
        );
      }

      //
      // 1. Solve Riemann problem (merge data)
      //
      solveRiemannProblemAtInterface(cellDescription1,cellDescription2,face);
    }

     cellDescription1.setNeighbourMergePerformed(face._faceIndex1,true);
     cellDescription2.setNeighbourMergePerformed(face._faceIndex2,true);
  }

  #ifdef USE_ITAC
  VT_end(mergeNeighboursHandle);
  #endif
}

std::string exahype::solvers::ADERDGSolver::riemannDataToString(
    const double* const Q,const double* const F,std::string suffix) const {
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  const int nodesPerFace = getBndFaceSize()/numberOfData;

  std::ostringstream stream;
  stream << std::endl;
  stream << "riemann states ("<<suffix<<"):" << std::endl;
  for (int i = 0; i<numberOfData; i++) {
    double minQ = std::numeric_limits<double>::infinity();
    double maxQ = -std::numeric_limits<double>::infinity();
    for(int n=0; n<nodesPerFace; ++n) {
      minQ = std::min(Q[n*numberOfData+i],minQ);
      maxQ = std::max(Q[n*numberOfData+i],maxQ);
    }
    stream << "Q"<<suffix<<"["<<i<<"] in ["<<std::setprecision(2)<<minQ<<","<<std::setprecision(2)<<maxQ<<"], ";
    stream << std::endl;
  }
  stream << "riemann fluxes ("<<suffix<<"):" << std::endl;
  for (int i = 0; i<_numberOfVariables; i++) {
    double minF = std::numeric_limits<double>::infinity();
    double maxF = -std::numeric_limits<double>::infinity();
    for(int n=0; n<nodesPerFace; ++n) {
      minF = std::min(F[n*_numberOfVariables+i],minF);
      maxF = std::max(F[n*_numberOfVariables+i],maxF);
    }
    stream << "F"<<suffix<<"["<<i<<"] in  ["<<std::setprecision(2)<<minF<<","<<std::setprecision(2)<<maxF<<"], ";
    stream << std::endl;
  }
  return stream.str();
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& cellDescription1,
    CellDescription& cellDescription2,
    Solver::InterfaceInfo& face) {
  CellDescription& pLeft  =
      (face._orientation1==1) ? cellDescription1 : cellDescription2;
  CellDescription& pRight =
      (face._orientation1==1) ? cellDescription2 : cellDescription1;

  assertion2(DataHeap::getInstance().isValidIndex(pLeft.getExtrapolatedPredictorIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pLeft.getFluctuationIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pRight.getExtrapolatedPredictorIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pRight.getFluctuationIndex()),pLeft.toString(),pRight.toString());

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace  = getBndFluxSize();

  double* QL = static_cast<double*>(pLeft.getExtrapolatedPredictor()) +  dataPerFace * face._faceIndexLeft;
  double* FL = static_cast<double*>(pLeft.getFluctuation()          ) +  dofsPerFace  * face._faceIndexLeft;

  double* QR = static_cast<double*>(pRight.getExtrapolatedPredictor()) + dataPerFace * face._faceIndexRight;
  double* FR = static_cast<double*>(pRight.getFluctuation()          ) + dofsPerFace  * face._faceIndexRight;

  // todo Time step must be interpolated in local time stepping case
  // both time step sizes are the same, so the min has no effect here.
  assertion3(std::isfinite(pLeft.getTimeStepSize()),pLeft.toString(),face._faceIndexLeft,face._direction);
  assertion3(std::isfinite(pRight.getTimeStepSize()),pRight.toString(),face._faceIndexRight,face._direction);
  assertion3(pLeft.getTimeStepSize()>=0.0,pLeft.toString(),face._faceIndexLeft,face._direction);
  assertion3(pRight.getTimeStepSize()>=0.0,pRight.toString(),face._faceIndexRight,face._direction);

  #ifdef Asserts
  std::string inputDataL = riemannDataToString(QL,FL,"L");
  std::string inputDataR = riemannDataToString(QR,FR,"R");
  #endif

  std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(pLeft,pRight);
  riemannSolver(
      FL,FR,QL,QR,
      std::get<0>(timeStepData),
      std::get<1>(timeStepData),
      pLeft.getSize(),
      face._direction, false, -1); // TODO(Dominic): Merge Riemann solver directly with the face integral and push the result on update
                                   // does not make sense to overwrite the flux when performing local time stepping; coarse grid flux must be constant, or not?

  #ifdef Asserts
  if ( _checkForNaNs ) { // assumes the solver is used as part of the hybrid solver
    std::string outputInformationL = riemannDataToString(QL,FL,"L");
    std::string outputInformationR = riemannDataToString(QR,FR,"R");

    const int nodesPerFace = dofsPerFace/_numberOfVariables;
    for (int i = 0; i<_numberOfVariables; i++) {
      for(int n=0; n<nodesPerFace; ++n) {
        assertion10(tarch::la::equals(pLeft.getTimeStepSize(),0.0) || (std::isfinite(FL[i]) && std::isfinite(FR[i])),
                   pLeft.toString(),pRight.toString(),face._direction,i,FL[i],FR[i],
                   inputDataL,inputDataR,outputInformationL,outputInformationR);
      }
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::mergeWithBoundaryData(CellDescription& cellDescription,const int faceIndex,const int direction,const int orientation) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  applyBoundaryConditions(cellDescription,faceIndex,direction,orientation);

  mergeWithAugmentationStatus(cellDescription,faceIndex,BoundaryStatus);
  mergeWithCommunicationStatus(cellDescription,faceIndex,BoundaryStatus);
  mergeWithRefinementStatus(cellDescription,faceIndex,BoundaryStatus);

  cellDescription.setNeighbourMergePerformed(faceIndex,true);
}

void exahype::solvers::ADERDGSolver::applyBoundaryConditions(CellDescription& cellDescription,const int faceIndex,const int direction,const int orientation) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()),cellDescription.toString());
  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("applyBoundaryConditions(...)","#boundaryConditions="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace = getBndFluxSize();
  const int gradientDataPerFace = getBndGradQSize();
  double* QIn = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) +  dataPerFace * faceIndex;
  double* gradQIn = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient()) +  gradientDataPerFace * faceIndex;
  double* FIn = static_cast<double*>(cellDescription.getFluctuation())           +  dofsPerFace * faceIndex;
  const double* luh = static_cast<double*>(cellDescription.getSolution());

  #ifdef Asserts
  std::string inputData = riemannDataToString(QIn,FIn,"In");
  #endif

  #ifdef Asserts
  assertion4(std::isfinite(cellDescription.getTimeStamp()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStamp());
  assertion4(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  assertion4(cellDescription.getTimeStepSize()>=0.0, cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  if ( _checkForNaNs ) {
    for(int i=0; i<dofsPerFace; ++i) {
      assertion6(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(FIn[i]),cellDescription.toString(),faceIndex,direction,i,FIn[i],inputData);
    }
  }
  #endif

  // TODO(Dominic): Hand in space-time volume data. Time integrate it afterwards

  std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(cellDescription,cellDescription);
  boundaryConditions(
      FIn,QIn, gradQIn,
      luh,
      cellDescription.getOffset() + 0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      std::get<0>(timeStepData),
      std::get<1>(timeStepData),
      direction,orientation);

  #ifdef Asserts
  assertion4(std::isfinite(cellDescription.getTimeStamp()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStamp());
  assertion4(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  assertion4(cellDescription.getTimeStepSize()>=0.0, cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  if ( _checkForNaNs ) {
    for(int i=0; i<dofsPerFace; ++i) {
      assertion6(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(FIn[i]),cellDescription.toString(),faceIndex,direction,i,FIn[i],inputData);
    }
  }
  #endif
}

#ifdef Parallel
const int exahype::solvers::ADERDGSolver::DataMessagesPerNeighbourCommunication    = 2;
const int exahype::solvers::ADERDGSolver::DataMessagesPerMasterWorkerCommunication = 2;

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void
exahype::solvers::ADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(cellDescription.getAugmentationStatus()); // TODO(Dominic): Add to docu: Might be merged multiple times!
    metadata.push_back(cellDescription.getCommunicationStatus());
    metadata.push_back(cellDescription.getRefinementStatus());
  } else {
    for (int i = 0; i < exahype::NeighbourCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const MetadataHeap::HeapEntries&             neighbourMetadata,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    const tarch::la::Vector<DIMENSIONS, int>&    posNeighbour,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre) {
  const int neighbourAugmentationStatus  = neighbourMetadata[exahype::NeighbourCommunicationMetadataAugmentationStatus  ];
  const int neighbourCommunicationStatus = neighbourMetadata[exahype::NeighbourCommunicationMetadataCommunicationStatus ];
  const int neighbourRefinementStatus    = neighbourMetadata[exahype::NeighbourCommunicationMetadataRefinementStatus    ];

  logDebug("mergeWithNeighbourMetadata(...)", "received neighbour metadata="<<neighbourAugmentationStatus<<","<<neighbourCommunicationStatus<<","<<neighbourRefinementStatus);

  mergeWithNeighbourMetadata(solverNumber,cellInfo,
      neighbourAugmentationStatus,neighbourCommunicationStatus,neighbourRefinementStatus,
      pos,posNeighbour,barycentre);
}

void exahype::solvers::ADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     solverNumber,
    Solver::CellInfo&                             cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  barycentre,
    const int                                     level) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    Solver::BoundaryFaceInfo face(src,dest);

    if ( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

      const int dofsPerFace = getBndFluxSize();
      const int dataPerFace = getBndFaceSize();

      const double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
      const double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation())           + dofsPerFace * face._faceIndex;

      waitUntilCompletedLastStep<CellDescription>(cellDescription,true,true);

      // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
      // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
      DataHeap::getInstance().sendData(
          lQhbnd, dataPerFace, toRank, barycentre, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().sendData(
          lFhbnd, dofsPerFace, toRank, barycentre, level,
          peano::heap::MessageType::NeighbourCommunication);
      // TODO(Dominic): If anarchic time stepping send the time step over too.
    }
  }
}

// TODO(Dominic): Add to docu: We only perform a Riemann solve if a Cell is involved.
void exahype::solvers::ADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const int                                    level) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      prefetchFaceData(cellDescription,face._faceIndex);

      // Send order: lQhbnd,lFhbnd
      // Receive order: lFhbnd,lQhbnd
      // TODO(Dominic): If anarchic time stepping, receive the time step too.
      const int dofsPerFace = getBndFluxSize();
      const int dataPerFace = getBndFaceSize();
      DataHeap::getInstance().receiveData(
          const_cast<double*>(_receivedFluctuations.data()),dofsPerFace, // TODO const-correct peano
          fromRank, barycentre, level,peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().receiveData(                              // TODO const-correct peano
          const_cast<double*>(_receivedExtrapolatedPredictor.data()),dataPerFace,
          fromRank, barycentre, level, peano::heap::MessageType::NeighbourCommunication);
      
      solveRiemannProblemAtInterface(
          cellDescription, face,
          _receivedExtrapolatedPredictor.data(),
          _receivedFluctuations.data(),
          fromRank);
    }

    cellDescription.setNeighbourMergePerformed(face._faceIndex,true);
  }
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& cellDescription,
    Solver::BoundaryFaceInfo& face,
    const double* const lQhbnd,
    const double* const lFhbnd,
    const int fromRank) {
  logDebug("solveRiemannProblemAtInterface(...)",
      "cell-description=" << cellDescription.toString());

  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace = getBndFluxSize();
  if ( face._orientation==0 ) {
    const double* const QL = lQhbnd;
    double* FL             = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    const double* const QR = static_cast<double*>( cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
    double* FR             = static_cast<double*>( cellDescription.getFluctuation())           + dofsPerFace * face._faceIndex;
    // TODO const-correct kernels

    #ifdef Asserts
    std::string inputDataL = riemannDataToString(QL,FL,"L");
    std::string inputDataR = riemannDataToString(QR,FR,"R");
    #endif

    riemannSolver(
        FL, FR, QL, QR,
	cellDescription.getTimeStamp(), cellDescription.getTimeStepSize(),
	cellDescription.getSize(), face._direction,false,face._faceIndex);
    #ifdef Asserts
    if ( _checkForNaNs ) {
      for (int ii = 0; ii<dofsPerFace; ii++) {
        assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
        assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
      }
    }
    #endif
  } else {
    const double* const QR = lQhbnd;
    const double* const QL = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
    double* FR = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    double* FL = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * face._faceIndex; // TODO const-correct kernels

    #ifdef Asserts
    std::string inputDataL = riemannDataToString(QL,FL,"L");
    std::string inputDataR = riemannDataToString(QR,FR,"R");
    #endif

    std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(cellDescription,cellDescription);
    riemannSolver(
        FL, FR, QL, QR,
        std::get<0>(timeStepData),
        std::get<1>(timeStepData),
	cellDescription.getSize(),
        face._direction,false,face._faceIndex);
    
    #ifdef Asserts
    if ( _checkForNaNs ) {
      for (int ii = 0; ii<dofsPerFace; ii++) {
        assertion10(std::isfinite(FL[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank,inputDataL,inputDataR);
        assertion10(std::isfinite(FR[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank,inputDataL,inputDataR);
      }
    }
    #endif
  }
}

void exahype::solvers::ADERDGSolver::dropNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries message = compileMessageForMaster();

  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("sendDataToMaster(...)","Sending data to master: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "to rank " << masterRank <<
        ", message size="<<message.size()
    );
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForMaster(const int capacity) const {
  const int messageSize = 2 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message;
  message.reserve(std::max(messageSize,capacity));

  message.push_back(_admissibleTimeStepSize);
  message.push_back(convertToDouble(_meshUpdateEvent));

  for (const auto observable : _globalObservables) {
    logDebug("sendDataToMaster(...)","Sending data to master: " << "entry=" << observable);
    message.push_back(observable);
  }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  assertion1(std::isfinite(message[0]),message[0]);
  return message;
}

/**
 * At the time of the merging,
 * the workers and the master have already performed
 * at local update of the next predictor time step size
 * and of the predictor time stamp.
 * We thus need to minimise over both quantities.
 */
void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 2 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithWorkerData(...)","Receive data from worker rank: " <<
             "data[0]=" << message[0] << "," <<
             "data[1]=" << message[1] << "," <<
             "from worker " << workerRank << "," <<
             "message size="<<message.size());
   }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  mergeWithWorkerData(message);

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithWorkerData(...)","Updated fields: " <<
             "_admissibleTimeStepSize=" << _admissibleTimeStepSize << "," <<
             "_meshUpdateEvent="        << Solver::toString(_meshUpdateEvent) );
  }
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerData(const DataHeap::HeapEntries& message) {
  int index=0; // post update
  _admissibleTimeStepSize = std::min( _admissibleTimeStepSize, message[index++] );
  _meshUpdateEvent       = mergeMeshUpdateEvents(_meshUpdateEvent,convertToMeshUpdateEvent(message[index++]));
  DataHeap::HeapEntries observablesFromWorker = DataHeap::HeapEntries(_numberOfGlobalObservables);
  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    observablesFromWorker[i] = message[index++];
  }
  mergeGlobalObservables(_nextGlobalObservables.data(), observablesFromWorker.data()); // !  master hasn't wrapped up time step yet
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForWorker(const int capacity) const {
  const auto messageSize = 5 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message;
  message.reserve(std::max(messageSize,capacity));

  message.push_back(_minTimeStamp);
  message.push_back(_minTimeStepSize);
  message.push_back(_estimatedTimeStepSize);
  message.push_back(convertToDouble(_meshUpdateEvent));
  message.push_back(_stabilityConditionWasViolated ? 1.0 : -1.0);

  for (const auto observable : _globalObservables) {
    message.push_back(observable);
  }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  return message;
}

void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries message = compileMessageForWorker();

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug(
        "sendDataToWorker(...)","Broadcasting time step data: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "data[2]=" << message[2] << "," <<
        "data[3]=" << message[3] << "," <<
        "data[4]=" << message[4]);
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(const DataHeap::HeapEntries& message) {
  int index=0;
  _minTimeStamp                  = message[index++];
  _minTimeStepSize               = message[index++];
  _estimatedTimeStepSize         = message[index++];
  _meshUpdateEvent               = convertToMeshUpdateEvent(message[index++]);
  _stabilityConditionWasViolated = (message[index++] > 0.0) ? true : false;

  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    _globalObservables[i] = message[index++];
  }
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 5 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  mergeWithMasterData(message);

  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug(
        "mergeWithMasterData(...)","Received message from master: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "data[2]=" << message[2] << "," <<
        "data[3]=" << message[3] << "," <<
        "data[4]=" << message[4]);
  }
}
#endif

std::string exahype::solvers::ADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::ADERDGSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:"                   << _identifier;
  out << "_type:"                         << Solver::toString(_type) << ",";
  out << "_numberOfVariables:"            << _numberOfVariables << ",";
  out << "_numberOfParameters:"           << _numberOfParameters << ",";
  out << "_nodesPerCoordinateAxis:"       << _nodesPerCoordinateAxis << ",";
  out << "_maximumMeshSize:"              << _maximumMeshSize << ",";
  out << "_timeStepping:"                 << Solver::toString(_timeStepping) << ",";
  out << "_unknownsPerFace:"              << getUnknownsPerFace() << ",";
  out << "_unknownsPerCellBoundary:"      << getUnknownsPerCellBoundary() << ",";
  out << "_unknownsPerCell:"              << getUnknownsPerCell() << ",";
  out << "_fluxUnknownsPerCell:"          << getFluxUnknownsPerCell() << ",";
  out << "_spaceTimeUnknownsPerCell:"     << getSpaceTimeUnknownsPerCell() << ",";
  out << "_spaceTimeFluxUnknownsPerCell:" << getSpaceTimeFluxUnknownsPerCell() << ",";
  out << "_previousMinTimeStamp:"         << _previousMinTimeStamp << ",";
  out << "_previousMinTimeStepSize:"      << _previousMinTimeStepSize << ",";
  out << "_minTimeStamp:"                 << _minTimeStamp << ",";
  out << "_minTimeStepSize:"              << _minTimeStepSize << ",";
  out << "_estimatedTimeStepSize:"        << _estimatedTimeStepSize << ",";
  out <<  ")";
}

exahype::solvers::ADERDGSolver::CompressionJob::CompressionJob(
  const ADERDGSolver& solver,
  CellDescription&    cellDescription,
  const bool          isSkeletonJob):
  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask,0,getTaskPriority(isSkeletonJob)),
  _solver(solver),
  _cellDescription(cellDescription),
  _isSkeletonJob(isSkeletonJob) {
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    NumberOfEnclaveJobs.fetch_add(1);
  }
}


bool exahype::solvers::ADERDGSolver::CompressionJob::run(bool runOnMasterThread) {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_sub(1);
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs.fetch_sub(1);
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
  return false;
}


void exahype::solvers::ADERDGSolver::compress( CellDescription& cellDescription, const bool isSkeletonCell ) const {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if ( SpawnCompressionAsBackgroundJob ) {
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);
      CompressionJob compressionJob( *this, cellDescription, isSkeletonCell );
      assertionMsg( false, "this call is invalid" );
/*
      if ( isSkeletonCell ) {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible  );
      } else {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::Background  );
      }
*/
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(CellDescription::Compressed);
    }
  }
}


void exahype::solvers::ADERDGSolver::uncompress(CellDescription& cellDescription) const {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    peano::datatraversal::TaskSet::finishToProcessBackgroundJobs();

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    madeDecision = cellDescription.getCompressionState() != CellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == CellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( CellDescription::CurrentlyProcessed );
    }
    lock.free();
  }
  #else
  bool uncompress = CompressionAccuracy>0.0
      && cellDescription.getCompressionState() == CellDescription::Compressed;
  #endif

  if (uncompress) {
    pullUnknownsFromByteStream(cellDescription);
    computeHierarchicalTransform(cellDescription,1.0);

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
      cellDescription.setCompressionState(CellDescription::Uncompressed);
    lock.free();
  }
}


void exahype::solvers::ADERDGSolver::determineUnknownAverages(
  CellDescription& cellDescription) const {
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdateAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationAveragesIndex()), cellDescription.toString() );

  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  double* solutionAverages              = static_cast<double*>(cellDescription.getSolutionAverages());
  double* previousSolutionAverage       = static_cast<double*>(cellDescription.getPreviousSolutionAverages());
  double* updateAverages                = static_cast<double*>(cellDescription.getUpdateAverages());
  double* extrapolatedPredictorAverages = static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages());
  double* fluctuationAverages           = static_cast<double*>(cellDescription.getFluctuationAverages());

  double* solution              = static_cast<double*>(cellDescription.getSolution());
  double* previousSolution      = static_cast<double*>(cellDescription.getPreviousSolution());
  double* update                = static_cast<double*>(cellDescription.getUpdate());
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* fluctuation           = static_cast<double*>(cellDescription.getFluctuation());

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      solutionAverages[variableNumber]        += solution        [idx_cellData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += previousSolution[idx_cellData(i,variableNumber)];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      updateAverages[variableNumber]          += update[idx_cellUnknowns(i,variableNumber)];
    }
  }
  for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
    solutionAverages[variableNumber]        = solutionAverages[variableNumber]        / (double) nodesPerCell;
    previousSolutionAverage[variableNumber] = previousSolutionAverage[variableNumber] / (double) nodesPerCell;
  }
  for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
    updateAverages[variableNumber]          = updateAverages[variableNumber]          / (double) nodesPerCell;
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
        extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] +=
            extrapolatedPredictor[idx_faceData(face,i,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
        fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] +=
            fluctuation[idx_faceUnknowns(face,i,variableNumber)];
      }
    }
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] =
          extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] / (double) nodesPerFace;
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] =
          fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)]       / (double) nodesPerFace;
    }
  }
}


void exahype::solvers::ADERDGSolver::computeHierarchicalTransform(
    CellDescription& cellDescription, double sign) const {
  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  double* solutionAverages              = static_cast<double*>(cellDescription.getSolutionAverages());
  double* previousSolutionAverage       = static_cast<double*>(cellDescription.getPreviousSolutionAverages());
  double* updateAverages                = static_cast<double*>(cellDescription.getUpdateAverages());
  double* extrapolatedPredictorAverages = static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages());
  double* fluctuationAverages           = static_cast<double*>(cellDescription.getFluctuationAverages());

  double* solution              = static_cast<double*>(cellDescription.getSolution());
  double* previousSolution      = static_cast<double*>(cellDescription.getPreviousSolution());
  double* update                = static_cast<double*>(cellDescription.getUpdate());
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* fluctuation           = static_cast<double*>(cellDescription.getFluctuation());

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      solution        [idx_cellData(i,variableNumber)] += sign * solutionAverages[variableNumber];
      previousSolution[idx_cellData(i,variableNumber)] += sign * previousSolutionAverage[variableNumber];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      update[idx_cellUnknowns(i,variableNumber)] += sign * updateAverages[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<DIMENSIONS_TIMES_TWO; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) {  // variables+parameters
        extrapolatedPredictor[idx_faceData(face,i,variableNumber)] +=
            sign * extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) {  // variables
        fluctuation[idx_faceUnknowns(face,i,variableNumber)] +=
            sign * fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)];
      }
    }
  }
}

void exahype::solvers::ADERDGSolver::putUnknownsIntoByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getUpdateCompressedIndex()==-1 );
  assertion( cellDescription.getExtrapolatedPredictorCompressedIndex()==-1 );
  assertion( cellDescription.getFluctuationCompressedIndex()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfUpdate;
  int compressionOfExtrapolatedPredictor;
  int compressionOfFluctuation;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdateIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuationIndex() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> bool { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getPreviousSolution()),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfUpdate = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getUpdate()),
      getUnknownsPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedPredictor = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getExtrapolatedPredictor()),
      getDataPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfFluctuation = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getFluctuation()),
      getUnknownsPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
    true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfUpdate);
  assertion(1<=compressionOfExtrapolatedPredictor);
  assertion(1<=compressionOfFluctuation);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfUpdate<=7);
  assertion(compressionOfExtrapolatedPredictor<=7);
  assertion(compressionOfFluctuation<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setPreviousSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getPreviousSolutionCompressedIndex()>=0 );
          cellDescription.setPreviousSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCell();
        tearApart(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), compressionOfPreviousSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionIndex(), true );
          cellDescription.setPreviousSolutionIndex(-1);
          cellDescription.setPreviousSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolutionCompressedIndex(CompressedDataHeap::getInstance().createData(0,0));
          assertion1( cellDescription.getSolutionCompressedIndex()>=0, cellDescription.getSolutionCompressedIndex() );
          cellDescription.setSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCell();

        tearApart(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), compressionOfSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getSolutionIndex(), true );
          cellDescription.setSolutionIndex(-1);
          cellDescription.setSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInUpdate(compressionOfUpdate);
      if (compressionOfUpdate<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setUpdateCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getUpdateCompressedIndex()>=0 );
          cellDescription.setUpdateCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getUpdateCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getUnknownsPerCell();
        tearApart(numberOfEntries, cellDescription.getUpdateIndex(), cellDescription.getUpdateCompressedIndex(), compressionOfUpdate);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getUpdateCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getUpdateIndex(), true );
          cellDescription.setUpdateIndex(-1);
          cellDescription.setUpdate(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInExtrapolatedPredictor(compressionOfExtrapolatedPredictor);
      if (compressionOfExtrapolatedPredictor<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setExtrapolatedPredictorCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getExtrapolatedPredictorCompressedIndex()>=0 );
          cellDescription.setExtrapolatedPredictorCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictorCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedPredictorIndex(), cellDescription.getExtrapolatedPredictorCompressedIndex(), compressionOfExtrapolatedPredictor);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorIndex(), true );
          cellDescription.setExtrapolatedPredictorIndex(-1);
          cellDescription.setExtrapolatedPredictor(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInFluctuation(compressionOfFluctuation);
      if (compressionOfFluctuation<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setFluctuationCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getFluctuationCompressedIndex()>=0 );
          cellDescription.setFluctuationCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getFluctuationCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getUnknownsPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getFluctuationIndex(), cellDescription.getFluctuationCompressedIndex(), compressionOfFluctuation);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getFluctuationCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getFluctuationIndex(), true );
          cellDescription.setFluctuationIndex(-1);
          cellDescription.setFluctuation(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}


void exahype::solvers::ADERDGSolver::pullUnknownsFromByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int dataPointsPerCell       = getDataPerCell();
  const int unknownsPerCellBoundary = getUnknownsPerCellBoundary();

  {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
      cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setSolutionIndex( DataHeap::getInstance().createData(         dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setUpdateIndex( DataHeap::getInstance().createData(           getUpdateSize(),   getUpdateSize() ) );

      cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      cellDescription.setFluctuationIndex( DataHeap::getInstance().createData(           unknownsPerCellBoundary, unknownsPerCellBoundary ) );
    lock.free();

    if (cellDescription.getPreviousSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getUpdateIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setUpdateIndex( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedPredictorIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData(unknownsPerCellBoundary ) );
      lock.free();
    }
    if (cellDescription.getFluctuationIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setFluctuationIndex( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      lock.free();
    }

    cellDescription.setPreviousSolution     ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionIndex()).data() ) );
    cellDescription.setSolution             ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionIndex()).data() ) );
    cellDescription.setExtrapolatedPredictor( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictorIndex()).data() ) );

  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ),
      cellDescription.getPreviousSolutionCompressedIndex()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ),
    cellDescription.getSolutionCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressedIndex() ),
    cellDescription.getUpdateCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressedIndex() ),
    cellDescription.getExtrapolatedPredictorCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressedIndex() ),
    cellDescription.getFluctuationCompressedIndex()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ), cellDescription.getPreviousSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressedIndex(), true );
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ), cellDescription.getSolutionIndex() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressedIndex(), true );
          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInUpdate()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getUpdateIndex() ), cellDescription.getUpdateIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressedIndex() ));
        const int numberOfEntries = getUnknownsPerCell();
        glueTogether(numberOfEntries, cellDescription.getUpdateIndex(), cellDescription.getUpdateCompressedIndex(), cellDescription.getBytesPerDoFInUpdate());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getUpdateCompressedIndex(), true );
          cellDescription.setUpdateCompressedIndex(-1);
          cellDescription.setUpdateCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInExtrapolatedPredictor()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorIndex() ), cellDescription.getExtrapolatedPredictorIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressedIndex() ));
        const int numberOfEntries = getDataPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedPredictorIndex(), cellDescription.getExtrapolatedPredictorCompressedIndex(), cellDescription.getBytesPerDoFInExtrapolatedPredictor());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorCompressedIndex(), true );
          cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
          cellDescription.setExtrapolatedPredictorCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInFluctuation()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuationIndex() ), cellDescription.getFluctuationIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressedIndex() ));
        const int numberOfEntries = getUnknownsPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getFluctuationIndex(), cellDescription.getFluctuationCompressedIndex(), cellDescription.getBytesPerDoFInFluctuation());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getFluctuationCompressedIndex(), true );
          cellDescription.setFluctuationCompressedIndex(-1);
          cellDescription.setFluctuationCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}

///////////////////////
// PROFILING
///////////////////////

exahype::solvers::Solver::CellProcessingTimes exahype::solvers::ADERDGSolver::measureCellProcessingTimes(const int numberOfRuns) {
  // Setup
  const int cellDescriptionsIndex = ADERDGSolver::Heap::getInstance().createData(0,1);
  FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionsIndex,0,1); // needs to be done

  Solver::CellInfo cellInfo(cellDescriptionsIndex);
  addNewCellDescription(
      0,cellInfo,CellDescription::Type::Leaf,
      getMaximumAdaptiveMeshLevel(), /* needs to be on the fine grid for the limiter cells */-1,
      getCoarsestMeshSize(),
      _domainOffset);

  CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[0];

  ensureNecessaryMemoryIsAllocated(cellDescription);

  adjustSolutionDuringMeshRefinementBody(cellDescription,true);
  updateTimeStepSize(0,cellInfo);
  const double dt = cellDescription.getTimeStepSize();

  // ADER-DG specific setup ( all Riemanns have been performed, cell is surrounded by other Cell type cells )
  cellDescription.setAugmentationStatus(0);
  cellDescription.setFacewiseAugmentationStatus(0);
  cellDescription.setCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);
  cellDescription.setFacewiseCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);

  // MEASUREMENTS
  CellProcessingTimes result;

  // measure ADERDG STP
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    int numberOfPicardIterations = std::numeric_limits<int>::max();
    for (int it=0; it<numberOfRuns; it++) {
      numberOfPicardIterations = predictionAndVolumeIntegralBody(cellDescription,cellDescription.getTimeStamp(),cellDescription.getTimeStepSize(),false,true,true);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._minTimePredictor = time_sec / numberOfRuns / std::abs(numberOfPicardIterations);
    if ( isLinear() ) {
      result._maxTimePredictor = result._minTimePredictor;
    } else {
      result._maxTimePredictor = result._minTimePredictor * getNodesPerCoordinateAxis(); // * (order+1)
    }
  }

  // measure ADER-DG cell update
  cellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
  cellDescription.setFacewiseRefinementStatus(_refineOrKeepOnFineGrid);
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      cellDescription.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> boundaryMarkers(0); // >= 0 indicates no remote/domain boundary
      updateBody(cellDescription,cellInfo,boundaryMarkers);

      swapSolutionAndPreviousSolution(cellDescription); // assumed to be very cheap
      rollbackToPreviousTimeStep(cellDescription);
      cellDescription.setTimeStepSize(dt);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure Riemann solve
  {
    const tarch::la::Vector<DIMENSIONS,int> pos1(0);
    tarch::la::Vector<DIMENSIONS,int> pos2(0); pos2[0]=1;
    Solver::InterfaceInfo face(pos1,pos2);
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solveRiemannProblemAtInterface(cellDescription,cellDescription,face);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGRiemann = time_sec / numberOfRuns;
  }

  // Clean up
  cellDescription.setType(CellDescription::Type::Erased);
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

  DataHeap::getInstance().deleteAllData();
  ADERDGSolver::Heap::getInstance().deleteAllData();
  FiniteVolumesSolver::Heap::getInstance().deleteAllData();

  return result;
}
