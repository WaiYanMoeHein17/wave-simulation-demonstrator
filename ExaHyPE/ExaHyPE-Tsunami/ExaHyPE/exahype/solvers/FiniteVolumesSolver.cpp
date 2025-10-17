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
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Güra
 **/

#include "exahype/solvers/FiniteVolumesSolver.h"

#include <iomanip>
#include <string>
#include <limits>
#include <algorithm> // copy_n
#include <chrono>    // profiling

#include "peano/utils/Loop.h"

#include "tarch/multicore/Lock.h"
#include "peano/datatraversal/TaskSet.h"

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "kernels/KernelUtils.h"

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "peano/heap/CompressedFloatingPointNumbers.h"

#include "tarch/multicore/Jobs.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"

namespace {
constexpr const char* tags[]{"solutionUpdate", "stableTimeStepSize"};
}  // namespace

#ifdef USE_ITAC
int exahype::solvers::FiniteVolumesSolver::adjustSolutionHandle     = 0;
int exahype::solvers::FiniteVolumesSolver::updateBodyHandle         = 0;
int exahype::solvers::FiniteVolumesSolver::updateBodyHandleSkeleton = 0;
int exahype::solvers::FiniteVolumesSolver::mergeNeighboursHandle    = 0;
#endif

tarch::logging::Log exahype::solvers::FiniteVolumesSolver::_log( "exahype::solvers::FiniteVolumesSolver");

void exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(const int cellDescriptionsIndex) {
  assertion(Heap::getInstance().isValidIndex(cellDescriptionsIndex));

  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    assertion(p.getType()==CellDescription::Type::Leaf);

    auto *solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

    FiniteVolumesSolver* fvSolver = nullptr;
    if (solver->getType()==Solver::Type::FiniteVolumes) {
      fvSolver = static_cast<FiniteVolumesSolver*>(solver);
    }
    else if (solver->getType()==Solver::Type::LimitingADERDG) {
      fvSolver =
          static_cast<LimitingADERDGSolver*>(solver)->getLimiter().get();
    }
    assertion(fvSolver!=nullptr);
    p.setType(CellDescription::Type::Erased);
    fvSolver->ensureNoUnnecessaryMemoryIsAllocated(p);
  }

  Heap::getInstance().getData(cellDescriptionsIndex).clear();
}

exahype::solvers::FiniteVolumesSolver::FiniteVolumesSolver(
    const std::string& identifier,
    const int numberOfVariables,
    const int numberOfParameters,
    const int numberOfGlobalObservables,
    const int basisSize,
    const int ghostLayerWidth,
    const double maximumMeshSize,
    const exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, exahype::solvers::Solver::Type::FiniteVolumes,
             numberOfVariables, numberOfParameters, numberOfGlobalObservables,
             basisSize,
             maximumMeshSize, 0,
             timeStepping, std::move(profiler)),
            _previousMinTimeStamp( std::numeric_limits<double>::infinity() ),
            _previousMinTimeStepSize( std::numeric_limits<double>::infinity() ),
            _minTimeStamp( std::numeric_limits<double>::infinity() ),
            _minTimeStepSize( std::numeric_limits<double>::infinity() ),
            _admissibleTimeStepSize( std::numeric_limits<double>::infinity() ),
            _ghostLayerWidth( ghostLayerWidth ),
            _meshUpdateEvent( MeshUpdateEvent::None ) {
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatch() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::FiniteVolumesSolver::getGhostLayerWidth() const {
  return _ghostLayerWidth;
}

int exahype::solvers::FiniteVolumesSolver::getGhostDataPerPatch() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis+2*_ghostLayerWidth, DIMENSIONS + 0) - getDataPerPatch();
}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatchFace() const {
  return _ghostLayerWidth*(_numberOfVariables+_numberOfParameters)*power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatchBoundary() const {
  return DIMENSIONS_TIMES_TWO *getDataPerPatchFace();
}

void exahype::solvers::FiniteVolumesSolver::resetMeshUpdateEvent() {
  _meshUpdateEvent = MeshUpdateEvent::None;
}

void exahype::solvers::FiniteVolumesSolver::updateMeshUpdateEvent(
    exahype::solvers::Solver::MeshUpdateEvent meshUpdateEvent) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _meshUpdateEvent = mergeMeshUpdateEvents(_meshUpdateEvent,meshUpdateEvent);
  lock.free();
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::FiniteVolumesSolver::getMeshUpdateEvent() const {
  return _meshUpdateEvent;
}

double exahype::solvers::FiniteVolumesSolver::getPreviousMinTimeStepSize() const {
  return _previousMinTimeStepSize;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}

void exahype::solvers::FiniteVolumesSolver::updateAdmissibleTimeStepSize(double value) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _admissibleTimeStepSize = std::min(_admissibleTimeStepSize, value);
  lock.free();
}

double exahype::solvers::FiniteVolumesSolver::getAdmissibleTimeStepSize() const {
  return _admissibleTimeStepSize;
}

void exahype::solvers::FiniteVolumesSolver::resetAdmissibleTimeStepSize() {
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
}

void exahype::solvers::FiniteVolumesSolver::initSolver(
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

  _previousMinTimeStepSize = 0.0;
  _previousMinTimeStamp = timeStamp;

  _minTimeStepSize = 0.0;
  _minTimeStamp = timeStamp;

  _meshUpdateEvent = MeshUpdateEvent::InitialRefinementRequested;

  // global observables
  _globalObservables.resize(_numberOfGlobalObservables);
  _nextGlobalObservables.resize(_numberOfGlobalObservables);
  resetGlobalObservables(_globalObservables.data()    );
  resetGlobalObservables(_nextGlobalObservables.data());
  
  // TODO(WORKAROUND)
  const int dataPerFace = getBndFaceSize();
  _invalidExtrapolatedSolution.resize(dataPerFace);
  std::fill_n(_invalidExtrapolatedSolution.data(),_invalidExtrapolatedSolution.size(),-1);

  init(cmdlineargs,parserView); // call user defined initalisation
}

bool exahype::solvers::FiniteVolumesSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  return false;
}

bool exahype::solvers::FiniteVolumesSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(CellDescription& cellDescription) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStamp);
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStepSize);
      cellDescription.setTimeStamp(_minTimeStamp);
      cellDescription.setTimeStepSize(_minTimeStepSize);
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::kickOffTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _meshUpdateEvent        = MeshUpdateEvent::None;
    _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
    resetGlobalObservables(_nextGlobalObservables.data());
  }

  // call user code
  beginTimeStep(_minTimeStamp,isFirstTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::FiniteVolumesSolver::wrapUpTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch,const bool isLastTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _previousMinTimeStepSize  = _minTimeStepSize;
    _previousMinTimeStamp     = _minTimeStamp;
  }
  _minTimeStamp += _minTimeStepSize;

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    _minTimeStepSize = _admissibleTimeStepSize;
  }

  if ( isLastTimeStepOfBatchOrNoBatch ) {
    std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),_globalObservables.begin());
    wrapUpGlobalObservables(_globalObservables.data());
  }

  // call user code
  endTimeStep(_minTimeStamp,isLastTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::FiniteVolumesSolver::updateTimeStepSize() {
  _minTimeStepSize        = _admissibleTimeStepSize;
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
}

void exahype::solvers::FiniteVolumesSolver::updateGlobalObservables() {
  std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),
            _globalObservables.begin());
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStep() {
  _minTimeStamp            = _previousMinTimeStamp;

  _previousMinTimeStamp    = std::numeric_limits<double>::infinity();
  _previousMinTimeStepSize = std::numeric_limits<double>::infinity();
}

bool exahype::solvers::FiniteVolumesSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) {
  bool result = cellDescriptionsIndex>=0;
  assertion1(!result || Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  return result;
}

int exahype::solvers::FiniteVolumesSolver::tryGetElement(
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

int exahype::solvers::FiniteVolumesSolver::computeWeight(const int cellDescriptionsIndex) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    return getCellDescriptions(cellDescriptionsIndex).size();
  }
  else return 0;
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  // Fine grid cell based uniform mesh refinement.
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
    fineGridCellElement==exahype::solvers::Solver::NotFound &&
    fineGridVerticesEnumerator.getLevel()==_coarsestMeshLevel
  ) {
    addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
               mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex,
               solverNumber);
    return true;
    // Fine grid cell based adaptive mesh refinement operations are not implemented.
  } else {
    return false;
  }
}

void exahype::solvers::FiniteVolumesSolver::addNewCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  CellInfo cellInfo =
      fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Leaf,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());

  const int fineGridCellElement = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[fineGridCellElement]; //TODO(Dominic): Multi-solvers: Might need to lock this?
  ensureNecessaryMemoryIsAllocated(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
    const int solverNumber,
    CellInfo& cellInfo,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Background job completion monitoring (must be initialised with true)
  newCellDescription.setHasCompletedLastStep(true);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  // newCellDescription.setHelperCellNeedsToStoreFaceData(false); // TODO(Dominic): Add to FV cell descr.

  newCellDescription.setNeighbourMergePerformed((signed char) 0/*implicit conversion*/);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Default data field indices
  newCellDescription.setSolutionIndex(-1);
  newCellDescription.setSolution(nullptr);
  newCellDescription.setSolutionAveragesIndex(-1);
  newCellDescription.setSolutionAverages(nullptr);
  newCellDescription.setSolutionCompressedIndex(-1);
  newCellDescription.setSolutionCompressed(nullptr);
  newCellDescription.setPreviousSolutionIndex(-1);
  newCellDescription.setPreviousSolution(nullptr);
  newCellDescription.setPreviousSolutionAveragesIndex(-1);
  newCellDescription.setPreviousSolutionAverages(nullptr);
  newCellDescription.setPreviousSolutionCompressedIndex(-1);
  newCellDescription.setPreviousSolutionCompressed(nullptr);
  newCellDescription.setExtrapolatedSolutionIndex(-1);
  newCellDescription.setExtrapolatedSolution(nullptr);
  newCellDescription.setExtrapolatedSolutionAveragesIndex(-1);
  newCellDescription.setExtrapolatedSolutionAverages(nullptr);
  newCellDescription.setExtrapolatedSolutionCompressedIndex(-1);
  newCellDescription.setExtrapolatedSolutionCompressed(nullptr);

  newCellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  cellInfo._FiniteVolumesCellDescriptions.push_back(newCellDescription);
  lock.free();
}


void exahype::solvers::FiniteVolumesSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())) {
    switch (cellDescription.getType()) {
      case CellDescription::Erased: {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionIndex()));

          if (cellDescription.getSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getSolutionAveragesIndex());

          if (cellDescription.getPreviousSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getPreviousSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAveragesIndex());

          if (cellDescription.getExtrapolatedSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getExtrapolatedSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionAveragesIndex());

          cellDescription.setSolutionIndex(-1);
          cellDescription.setSolution(nullptr);
          cellDescription.setPreviousSolutionIndex(-1);
          cellDescription.setPreviousSolution(nullptr);
          cellDescription.setExtrapolatedSolutionIndex(-1);
          cellDescription.setExtrapolatedSolution(nullptr);

          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);
          cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
          cellDescription.setExtrapolatedSolutionCompressed(nullptr);

          cellDescription.setSolutionAveragesIndex(-1);
          cellDescription.setSolutionAverages(nullptr);
          cellDescription.setPreviousSolutionAveragesIndex(-1);
          cellDescription.setPreviousSolutionAverages(nullptr);
          cellDescription.setExtrapolatedSolutionAveragesIndex(-1);
          cellDescription.setExtrapolatedSolutionAverages(nullptr);
        lock.free();
      } break;
      case CellDescription::Type::Leaf:
        // do nothing
        break;
      default:
        assertionMsg(false,"No other cell description types are supported at the moment!");
        break;
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::checkDataHeapIndex(const CellDescription& cellDescription, const int arrayIndex,const std::string arrayName) {
  assertion1(DataHeap::getInstance().isValidIndex(arrayIndex),cellDescription.toString());
  if ( arrayIndex < 0 ) {
    logError("checkDataHeapIndex(...)","The data heap array 'cellDescription."<<arrayName<<"' could not be allocated! Potential reason: Not enough memory available." <<
             " CellDescription="<<cellDescription.toString());
    std::abort();
  }
}

void exahype::solvers::FiniteVolumesSolver::ensureNecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  switch (cellDescription.getType()) {
    case CellDescription::Type::Leaf:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
        // Allocate volume data
        const int patchSize = getDataPerPatch()+getGhostDataPerPatch();

        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolutionIndex(        DataHeap::getInstance().createData( patchSize, patchSize ));
          cellDescription.setPreviousSolutionIndex(DataHeap::getInstance().createData( patchSize, patchSize ));
          checkDataHeapIndex(cellDescription,cellDescription.getSolutionIndex(),"getSolutionIndex()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionIndex(),"getPreviousSolutionIndex()");
          cellDescription.setSolution        ( getDataHeapEntries(cellDescription.getSolutionIndex()).data() ) ;
          cellDescription.setPreviousSolution( getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).data() ) ;
          // std::fill_n(static_cast<double*>(cellDescription.getPreviousSolution()),patchSize,std::numeric_limits<double>::quiet_NaN());
          // std::fill_n(static_cast<double*>(cellDescription.getSolution()),patchSize,std::numeric_limits<double>::quiet_NaN());
          // Zero out the solution and previous solution arrays. For our MUSCL-Hancock implementation which
          // does not take the corner neighbours into account e.g., it is important that the values in
          // the corner cells of the first ghost layer are set to zero.
          std::fill_n( static_cast<double*>(cellDescription.getSolution()),         patchSize, 0.0 );
          std::fill_n( static_cast<double*>(cellDescription.getPreviousSolution()), patchSize, 0.0 );

          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);

          const int dataPerSubcell = getNumberOfVariables()+getNumberOfParameters();
          cellDescription.setSolutionAveragesIndex        (DataHeap::getInstance().createData( dataPerSubcell, dataPerSubcell ) );
          cellDescription.setPreviousSolutionAveragesIndex(DataHeap::getInstance().createData( dataPerSubcell, dataPerSubcell ) );
          checkDataHeapIndex(cellDescription,cellDescription.getSolutionAveragesIndex(),"getSolutionAveragesIndex()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionAveragesIndex(),"getPreviousSolutionAveragesIndex()");
          cellDescription.setSolutionAverages        ( getDataHeapEntries(cellDescription.getSolutionAveragesIndex()).data() ) ;
          cellDescription.setPreviousSolutionAverages( getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex()).data() ) ;
          std::fill_n(static_cast<double*>(cellDescription.getPreviousSolutionAverages()),dataPerSubcell,std::numeric_limits<double>::quiet_NaN());
          std::fill_n(static_cast<double*>(cellDescription.getSolutionAverages()),dataPerSubcell,std::numeric_limits<double>::quiet_NaN());

          // Allocate boundary data
          const int patchBoundarySize = getDataPerPatchBoundary();
          cellDescription.setExtrapolatedSolutionIndex(DataHeap::getInstance().createData( patchBoundarySize, patchBoundarySize ));
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolutionIndex(),"getExtrapolatedSolutionIndex()");
          cellDescription.setExtrapolatedSolution( getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).data() ) ;
          std::fill_n( static_cast<double*>(cellDescription.getExtrapolatedSolution()), patchBoundarySize, 0.0 );

          cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
          cellDescription.setExtrapolatedSolutionCompressed(nullptr);

          cellDescription.setExtrapolatedSolutionAveragesIndex( DataHeap::getInstance().createData(dataPerSubcell*2*DIMENSIONS, dataPerSubcell*2*DIMENSIONS ) );
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolutionAveragesIndex(),"getExtrapolatedSolutionAveragesIndex()");
          cellDescription.setExtrapolatedSolutionAverages( getDataHeapEntries(cellDescription.getExtrapolatedSolutionAveragesIndex()).data() ) ;
          std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedSolutionAverages()),dataPerSubcell*2*DIMENSIONS,std::numeric_limits<double>::quiet_NaN());
        lock.free();
      }
      break;
    case CellDescription::Erased:
      // do nothing
    default:
      assertionMsg(false,"No other cell description types are supported at the moment!");
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
     const int solverNumber) {
  // do nothing
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::FiniteVolumesSolver::eraseOrRefineAdjacentVertices(
    const int cellDescriptionsIndex,
    const int solverNumber,
    const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const bool checkThoroughly,
    bool& checkSuccessful) const {
  checkSuccessful = true;
  if ( level < _coarsestMeshLevel ) {
    return RefinementControl::Refine;
  } else {
    return RefinementControl::Erase;
  }
}

void exahype::solvers::FiniteVolumesSolver::finaliseStateUpdates(
      const int solverNumber,
      CellInfo& cellInfo) {
  // do nothing
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////

double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(CellDescription& cellDescription,const bool isFirstTimeStepOfBatch) {
  assertion1(cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Type::Leaf,cellDescription.toString());
  // n-1
  if (isFirstTimeStepOfBatch) {
    cellDescription.setPreviousTimeStamp(cellDescription.getTimeStamp());
    cellDescription.setPreviousTimeStepSize(cellDescription.getTimeStepSize());
  }
  // n
  cellDescription.setTimeStamp(cellDescription.getTimeStamp()+cellDescription.getTimeStepSize());
  if ( getTimeStepping() != TimeStepping::GlobalFixed ) {
    double* solution = static_cast<double*>(cellDescription.getSolution());
    double admissibleTimeStepSize = stableTimeStepSize(solution, cellDescription.getSize());
    assertion1(std::isfinite(admissibleTimeStepSize),cellDescription.toString());
    cellDescription.setTimeStepSize(admissibleTimeStepSize);
    return admissibleTimeStepSize;
  } else {
    return cellDescription.getTimeStepSize();
  }
}

void exahype::solvers::FiniteVolumesSolver::updateTimeStepSize(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      double* solution = static_cast<double*>(cellDescription.getSolution());
      double admissibleTimeStepSize = stableTimeStepSize(solution, cellDescription.getSize());
      assertion1(std::isfinite(admissibleTimeStepSize),cellDescription.toString());
      cellDescription.setTimeStepSize(admissibleTimeStepSize);
      updateAdmissibleTimeStepSize(admissibleTimeStepSize);
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::updateGlobalObservables(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];

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
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStep(CellDescription& cellDescription) const {
  cellDescription.setTimeStamp(cellDescription.getPreviousTimeStamp());
  cellDescription.setTimeStepSize(cellDescription.getPreviousTimeStepSize());

  cellDescription.setPreviousTimeStamp(std::numeric_limits<double>::infinity());
  cellDescription.setPreviousTimeStepSize(std::numeric_limits<double>::infinity());
}

void exahype::solvers::FiniteVolumesSolver::adjustSolutionDuringMeshRefinementBody(
    CellDescription& cellDescription,
    const bool isInitialMeshRefinement) {
  #ifdef USE_ITAC
  VT_begin(adjustSolutionHandle);
  #endif

  assertion(cellDescription.getType()==CellDescription::Type::Leaf);

  adjustSolution(cellDescription);

  cellDescription.setHasCompletedLastStep(true);

  #ifdef USE_ITAC
  VT_end(adjustSolutionHandle);
  #endif
}

void exahype::solvers::FiniteVolumesSolver::adjustSolution(CellDescription& cellDescription) {
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
  for (int i=0; i<getDataPerPatch()+getGhostDataPerPatch(); i++) {
    assertion3(std::isfinite(solution[i]),cellDescription.toString(),"setInitialConditions(...)",i);
  }
  #endif
}

void exahype::solvers::FiniteVolumesSolver::reduce(
    const CellDescription& cellDescription,
    const UpdateResult&    result) {
  updateAdmissibleTimeStepSize(result._timeStepSize);

  if ( _numberOfGlobalObservables > 0 ) {
    assert(cellDescription.getType()==CellDescription::Type::Leaf);
    const double* const luh = static_cast<double*>(cellDescription.getSolution());
    const auto& cellSize    = cellDescription.getSize();
    const auto cellCentre   = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
    const auto t             = cellDescription.getTimeStamp();
    const auto dt            = cellDescription.getTimeStepSize();
    updateGlobalObservables(_nextGlobalObservables.data(),luh,cellCentre,cellSize,t,dt);
  }
}

void exahype::solvers::FiniteVolumesSolver::updateBody(
    CellDescription&                                   cellDescription,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         uncompressBefore) {
  #ifdef USE_ITAC
  if ( isAtRemoteBoundary ) {
    VT_begin(updateBodyHandleSkeleton);
  } else {
    VT_begin(updateBodyHandle);
  }
  #endif

  if ( uncompressBefore ) { uncompress(cellDescription); }

  updateSolution(cellDescription,cellInfo._cellDescriptionsIndex,boundaryMarkers,isFirstTimeStepOfBatch);
  UpdateResult result;
  result._timeStepSize = startNewTimeStep(cellDescription,isFirstTimeStepOfBatch);

  reduce(cellDescription,result);

  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
  compress(cellDescription,isAtRemoteBoundary);

  cellDescription.setHasCompletedLastStep(true); // last step of the FV update

  #ifdef USE_ITAC
  if ( isAtRemoteBoundary ) {
    VT_end(updateBodyHandleSkeleton);
  } else {
    VT_end(updateBodyHandle);
  }
  #endif
}

void exahype::solvers::FiniteVolumesSolver::fusedTimeStepOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);

  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if (
      (element != NotFound) &&
      (SpawnUpdateAsBackgroundJob || (SpawnPredictionAsBackgroundJob && !isLastTimeStepOfBatch))
  ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    const bool isSkeletonCell = isAtRemoteBoundary;
    peano::datatraversal::TaskSet( new FusedTimeStepJob(
        *this, cellDescription, cellInfo,
        isFirstTimeStepOfBatch, isLastTimeStepOfBatch,
        boundaryMarkers,isSkeletonCell ) );
  }
  else if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);
    updateBody(
        cellDescription,cellInfo,
        isFirstTimeStepOfBatch,isLastTimeStepOfBatch,boundaryMarkers,false/*uncompressBefore*/);
  }
}

void exahype::solvers::FiniteVolumesSolver::updateOrRestrict(
      const int                                          solverNumber,
      CellInfo&                                          cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers){
  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);

  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( (element!=NotFound) && SpawnUpdateAsBackgroundJob) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    peano::datatraversal::TaskSet(
        new UpdateJob(*this,cellDescription,cellInfo,boundaryMarkers) );
  }
  else if ( element!=NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    updateBody(
        cellDescription,cellInfo,
        true,true,boundaryMarkers,true/*uncompressBefore*/);
  }
}

void exahype::solvers::FiniteVolumesSolver::compress(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    compress(cellDescription,isAtRemoteBoundary);
  }
}

void exahype::solvers::FiniteVolumesSolver::adjustSolutionDuringMeshRefinement(
    const int solverNumber,CellInfo& cellInfo) {
  const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);

    if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
      peano::datatraversal::TaskSet( new AdjustSolutionDuringMeshRefinementJob(*this,cellDescription,isInitialMeshRefinement) );
    } else {
      adjustSolutionDuringMeshRefinementBody(cellDescription,isInitialMeshRefinement);
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    CellDescription&                                   cellDescription,
    const int                                          cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         backupPreviousSolution) {
  // boundary conditions
  double* solution       = static_cast<double*>(cellDescription.getSolution());
  for (int direction=0; direction<DIMENSIONS; direction++) {
    for (int orientation=0; orientation<2; orientation++) {
      const int faceIndex=2*direction+orientation;
      if ( boundaryMarkers[faceIndex]==mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex ) {
        tarch::la::Vector<DIMENSIONS,int> posCell(1);
        tarch::la::Vector<DIMENSIONS,int> posBoundary = posCell;
        posBoundary[direction] += orientation > 0 ? 1 : -1;

        boundaryConditions(
            solution,
            cellDescription.getOffset()+0.5*cellDescription.getSize(),
            cellDescription.getSize(),
            cellDescription.getTimeStamp(),
            cellDescription.getTimeStepSize(),
            posCell,posBoundary);

        cellDescription.setNeighbourMergePerformed(faceIndex,true);
      }
    }
  }

  if ( !tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) && !ProfileUpdate ) {
    logError("updateSolution(...)","Not all ghost layers were copied to cell="<<cellDescription.toString());
    std::terminate();
  }
  assertion1( tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) || ProfileUpdate,cellDescription.toString());
  cellDescription.setNeighbourMergePerformed(static_cast<signed char>(false));

  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
    static int counter = 0;
    static double myTimeStamp = 0;
    if ( !tarch::la::equals(myTimeStamp,_minTimeStamp,1e-9) ) {
      logInfo("mergeNeighboursData(...)","#updateSolution="<<counter);
      myTimeStamp = _minTimeStamp;
      counter=0;
    }
    counter++;
  #endif

  double* solutionBackup = static_cast<double*>(cellDescription.getPreviousSolution());
  if ( backupPreviousSolution ) {
    std::copy(solution,solution+getDataPerPatch()+getGhostDataPerPatch(),solutionBackup); // Copy (current solution) in old solution field.
  }

  assertion3(cellDescription.getTimeStepSize() > 0 && std::isfinite(cellDescription.getTimeStamp()) && std::isfinite(cellDescription.getTimeStepSize()), cellDescription.getTimeStamp(),cellDescription.getTimeStepSize(),cellDescription.toString());
  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[pre]");

  double admissibleTimeStepSize=std::numeric_limits<double>::infinity();

  // update
  solutionUpdate(solution,
		 cellDescription.getOffset() + 0.5*cellDescription.getSize(),
		 cellDescription.getSize(),
		 cellDescription.getTimeStamp(),
		 cellDescription.getTimeStepSize(),
		 admissibleTimeStepSize);

  if ( !tarch::la::equals(cellDescription.getTimeStepSize(), 0.0) && tarch::la::smaller(admissibleTimeStepSize,cellDescription.getTimeStepSize()) ) {
    logWarning("updateSolution(...)","Finite volumes solver time step size harmed CFL condition. dt="<<
               cellDescription.getTimeStepSize()<<", dt_adm=" << admissibleTimeStepSize << ". cell=" <<cellDescription.toString());
  }

  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),
      cellDescription.getTimeStepSize());

  // only for profiling
  if ( Solver::ProfileUpdate ) { swapSolutionAndPreviousSolution(cellDescription); }

  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[post]");
}

void exahype::solvers::FiniteVolumesSolver::swapSolutionAndPreviousSolution(
    CellDescription& cellDescription) const {
  // Simply swap the heap indices
  const int previousSolutionIndex = cellDescription.getPreviousSolutionIndex();
  void* previousSolution          = cellDescription.getPreviousSolution(); // pointer
  cellDescription.setPreviousSolutionIndex(cellDescription.getSolutionIndex());
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolutionIndex(previousSolutionIndex);
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::FiniteVolumesSolver::rollbackSolutionGlobally(const int solverNumber,CellInfo& cellInfo) const {
  // do nothing
  logError("rollbackSolutionGlobally(...)","Not implemented");
  std::abort();
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::mergeNeighboursData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo1,
    Solver::CellInfo&                         cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  #ifdef USE_ITAC
  VT_begin(mergeNeighboursHandle);
  #endif

  const int element1 = cellInfo1.indexOfFiniteVolumesCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    CellDescription& cellDescription1 = cellInfo1._FiniteVolumesCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._FiniteVolumesCellDescriptions[element2];

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

    waitUntilCompletedLastStep<CellDescription>(cellDescription1,false,false);
    waitUntilCompletedLastStep<CellDescription>(cellDescription2,false,false);

    assertion(cellDescription1.getType()==CellDescription::Type::Leaf && cellDescription2.getType()==CellDescription::Type::Leaf);

    assertion2(cellDescription1.getTimeStamp()<std::numeric_limits<double>::infinity(),cellDescription1.toString(),cellInfo1._cellDescriptionsIndex);
    assertion2(cellDescription1.getTimeStepSize()<std::numeric_limits<double>::infinity(),cellDescription1.toString(),cellInfo1._cellDescriptionsIndex);
    assertion2(cellDescription2.getTimeStamp()<std::numeric_limits<double>::infinity(),cellDescription2.toString(),cellInfo2._cellDescriptionsIndex);
    assertion2(cellDescription2.getTimeStepSize()<std::numeric_limits<double>::infinity(),cellDescription2.toString(),cellInfo2._cellDescriptionsIndex);

    if ( CompressionAccuracy > 0.0 ) {
     peano::datatraversal::TaskSet uncompression(
         [&] () -> bool {
       uncompress(cellDescription1);
       return false;
     },
     [&] () -> bool {
       uncompress(cellDescription2);
       return false;
     },
     peano::datatraversal::TaskSet::TaskType::Background,
     peano::datatraversal::TaskSet::TaskType::Background,
     true
     );
    }

    double* solution1 = static_cast<double*>(cellDescription1.getSolution());
    double* solution2 = static_cast<double*>(cellDescription2.getSolution());

    ghostLayerFilling(solution1,solution2,pos2-pos1);
    ghostLayerFilling(solution2,solution1,pos1-pos2);

    Solver::InterfaceInfo face(pos1,pos2);
    cellDescription1.setNeighbourMergePerformed(face._faceIndex1,true);
    cellDescription2.setNeighbourMergePerformed(face._faceIndex2,true);
  }

  #ifdef USE_ITAC
  VT_end(mergeNeighboursHandle);
  #endif
}

#ifdef Parallel
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerForkOrJoinCommunication   = 2;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerMasterWorkerCommunication = 0;

void exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    for (auto& cellDescription : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if ( !cellDescription.getHasCompletedLastStep() ) {
        peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
      }
      int numberOfBackgroundJobsToProcess = 1;
      while ( !cellDescription.getHasCompletedLastStep() ) {
        tarch::multicore::jobs::processBackgroundJobs(numberOfBackgroundJobsToProcess);
        numberOfBackgroundJobsToProcess++;
      }
    }
    Heap::getInstance().sendData(cellDescriptionsIndex,toRank,x,level,messageType);
  } else {
    sendEmptyCellDescriptions(toRank,messageType,x,level);
  }
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::FiniteVolumesSolver::ensureOnlyNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  auto* solver = RegisteredSolvers[cellDescription.getSolverNumber()];
  switch (solver->getType()) {
  case exahype::solvers::Solver::Type::ADERDG:
    assertionMsg(false,"Solver type not supported!");
    break;
  case exahype::solvers::Solver::Type::LimitingADERDG:
    static_cast<LimitingADERDGSolver*>(solver)->
        getLimiter()->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<LimitingADERDGSolver*>(solver)->
        getLimiter()->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::FiniteVolumes:
    static_cast<FiniteVolumesSolver*>(solver)->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<FiniteVolumesSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  }
}

void exahype::solvers::FiniteVolumesSolver::receiveCellDescriptions(
    const int                                     fromRank,
    exahype::Cell&                                localCell,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(
      localCell.getCellDescriptionsIndex(),fromRank,x,level,messageType);

  logDebug("mergeCellDescriptionsWithRemoteData(...)","received " <<
          Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).size() <<
          " cell descriptions for cell (centre="<< x.toString() << "level="<< level << ")");

  for (auto& cellDescription : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
    resetIndicesAndFlagsOfReceivedCellDescription(
        cellDescription,mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
  }
}

void exahype::solvers::FiniteVolumesSolver::resetIndicesAndFlagsOfReceivedCellDescription(
    CellDescription& cellDescription,const int parentIndex) {
  cellDescription.setParentIndex(parentIndex);

  // Default field data indices
  cellDescription.setSolutionIndex(-1);
  cellDescription.setSolution(nullptr);
  cellDescription.setPreviousSolutionIndex(-1);
  cellDescription.setPreviousSolution(nullptr);
  cellDescription.setExtrapolatedSolutionIndex(-1);
  cellDescription.setExtrapolatedSolution(nullptr);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
  cellDescription.setExtrapolatedSolutionCompressed(nullptr);
  cellDescription.setSolutionCompressedIndex(-1);
  cellDescription.setSolutionCompressed(nullptr);
  cellDescription.setPreviousSolutionCompressedIndex(-1);
  cellDescription.setPreviousSolutionCompressed(nullptr);

  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setExtrapolatedSolutionAveragesIndex(-1);
  cellDescription.setExtrapolatedSolutionAverages(nullptr);

  cellDescription.setBytesPerDoFInPreviousSolution(-1);
  cellDescription.setBytesPerDoFInSolution(-1);
  cellDescription.setBytesPerDoFInExtrapolatedSolution(-1);

  // background jobs
  cellDescription.setHasCompletedLastStep(true);
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::FiniteVolumesSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

////////////////////////////////////
// MASTER <=> WORKER
////////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
      element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << cellDescription.getSolverNumber() << " sent to rank "<<toRank<<
      ", cell: "<< x << ", level: " << level);

  assertion(cellDescription.getType()==CellDescription::Type::Leaf);
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),
      cellDescriptionsIndex,cellDescription.toString());
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()),
      cellDescriptionsIndex,cellDescription.toString());
  DataHeap::getInstance().sendData(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level, messageType);
  DataHeap::getInstance().sendData(
      static_cast<double*>(cellDescription.getPreviousSolution()),
      getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level, messageType);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  auto& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);
  assertion(cellDescription.getType()==CellDescription::Type::Leaf);

  logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
           ", cell: "<< x << ", level: " << level);

  // allocate memory
  ensureNecessaryMemoryIsAllocated(cellDescription);
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

  getDataHeapEntries(cellDescription.getSolutionIndex()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getSolutionIndex(),
      fromRank, x, level, messageType);
  getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getPreviousSolutionIndex(),
      fromRank, x, level, messageType );
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerIfProlongating(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::receiveDataFromMasterIfProlongating(
    const int masterRank,
    const int receivedCellDescriptionsIndex,
    const int receivedElement,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  // do nothing
}

bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
 // do nothing
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  // do nothing
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void
exahype::solvers::FiniteVolumesSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription =
        getCellDescription(cellDescriptionsIndex,element);
    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(exahype::InvalidMetadataEntry);
    metadata.push_back(exahype::InvalidMetadataEntry);
    metadata.push_back(exahype::InvalidMetadataEntry);
  } else {
    for (int i = 0; i < exahype::NeighbourCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     solverNumber,
    Solver::CellInfo&                             cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  barycentre,
    const int                                     level) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    Solver::BoundaryFaceInfo face(src,dest);

    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    waitUntilCompletedLastStep<CellDescription>(cellDescription,true,true);

    const int dataPerFace = getDataPerPatchFace();
    double* luhbnd    = static_cast<double*>(cellDescription.getExtrapolatedSolution()) + (face._faceIndex * dataPerFace);
    const double* luh = static_cast<double*>(cellDescription.getSolution());
    boundaryLayerExtraction(luhbnd,luh,dest-src);

    logDebug("sendDataToNeighbour(...)","send "<<DataMessagesPerNeighbourCommunication<<" arrays to rank=" <<toRank << ",cell="<<cellDescription.getOffset()<<",barycentre="<<barycentre<<",level="<<level);

    DataHeap::getInstance().sendData(
        luhbnd, dataPerFace, toRank, barycentre, level,
        peano::heap::MessageType::NeighbourCommunication);
  }
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  barycentre,
    const int                                     level) const {
  // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
  // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
  // TODO(WORKAROUND)
  #if defined(UsePeanosSymmetricBoundaryExchanger)
  DataHeap::getInstance().sendData(
      _invalidExtrapolatedSolution, toRank, barycentre, level,
      peano::heap::MessageType::NeighbourCommunication);
  #else
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        exahype::EmptyDataHeapMessage, toRank, barycentre, level,
        peano::heap::MessageType::NeighbourCommunication);
  #endif
}

void exahype::solvers::FiniteVolumesSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const int                                    level) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];

    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    logDebug("mergeWithNeighbourData(...)", "receive "<<DataMessagesPerNeighbourCommunication<<" arrays from rank="<<fromRank<<",x="<<barycentre<<",level="<<level);

    // TODO(Dominic): If anarchic time stepping, receive the time step too.
    //
    // Copy the received boundary layer into a ghost layer of the solution.
    // TODO(Dominic): Pipe it directly through the Riemann solver if
    // we only use the Godunov method and not higher-order FVM methods.
    // For methods that are higher order in time, e.g., MUSCL-Hancock, we usually need
    // corner neighbours. This is why we currently adapt a GATHER-UPDATE algorithm
    // instead of a SOLVE RIEMANN PROBLEM AT BOUNDARY-UPDATE INTERIOR scheme.
    const int dataPerFace = getDataPerPatchFace();
    double* luhbnd = static_cast<double*>(cellDescription.getExtrapolatedSolution()) + (face._faceIndex * dataPerFace);
 
    DataHeap::getInstance().receiveData(
        luhbnd, dataPerFace, fromRank, barycentre, level,peano::heap::MessageType::NeighbourCommunication);

    double* luh = static_cast<double*>(cellDescription.getSolution());
    ghostLayerFillingAtBoundary(luh,luhbnd,src-dest);

    cellDescription.setNeighbourMergePerformed(face._faceIndex,true);
  }
}

void exahype::solvers::FiniteVolumesSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  barycentre,
    const int                                     level) const {
  logDebug("dropNeighbourData(...)", "drop "<<DataMessagesPerNeighbourCommunication<<" arrays from rank="<<fromRank<<",x="<<barycentre<<",level="<<level);

  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, barycentre, level,
        peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const auto messageSize = 1 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message;
  message.reserve(messageSize);
  message.push_back(_admissibleTimeStepSize);
  for (const auto observable : _globalObservables) {
    message.push_back(observable);
  }

  assertion1(std::isfinite(message[0]),message[0]);
  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("sendDataToMaster(...)","Sending time step data: " << "data[0]=" << message[0]);
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 1 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(std::isfinite(message[0]),message[0]);

  int index=0;
  _admissibleTimeStepSize = std::min( _admissibleTimeStepSize, message[index++] );
  DataHeap::HeapEntries observablesFromWorker = DataHeap::HeapEntries(_numberOfGlobalObservables);
  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    observablesFromWorker[i] = message[index++];
  }
  mergeGlobalObservables(_nextGlobalObservables.data(), observablesFromWorker.data()); // !  master hasn't wrapped up time step yet

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithWorkerData(...)","Received data: " <<
             "data[0]=" << message[0]);
    logDebug("mergeWithWorkerData(...)","Updated fields: " <<
             "_admissibleTimeStepSize=" << _admissibleTimeStepSize);
  }
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  std::vector<double> message(0,2);
  message.push_back(_minTimeStamp);
  message.push_back(_minTimeStepSize);
  for (const auto observable : _globalObservables) {
    message.push_back(observable);
  }

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("sendDataToWorker(...)","Broadcasting time step data: " <<
        " data[0]=" << message[0] <<
        ",data[1]=" << message[1]);
    logDebug("sendDataWorker(...)","_minNextTimeStepSize="<<_admissibleTimeStepSize);
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 5 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithMasterData(...)","Received message from master: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1]);
  }

  int index = 0;
  _minTimeStamp    = message[index++];
  _minTimeStepSize = message[index++];
  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    _globalObservables[i] = message[index++];
  }
}
#endif

void exahype::solvers::FiniteVolumesSolver::validateNoNansInFiniteVolumesSolution(
    CellDescription& cellDescription,const int cellDescriptionsIndex,const char* methodTrace)  const {
  #if defined(Asserts)
  double* solution = static_cast<double*>(cellDescription.getSolution());

  dfor(i,_nodesPerCoordinateAxis+_ghostLayerWidth) {
    if (tarch::la::allSmaller(i,_nodesPerCoordinateAxis+_ghostLayerWidth)
    && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
      for (int unknown=0; unknown < _numberOfVariables; unknown++) {
        int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
        // cellDescription.getTimeStepSize()==0.0 is an initial condition
        assertion7(std::isfinite(solution[iScalar]),
                   cellDescription.toString(),cellDescriptionsIndex,solution[iScalar],i.toString(),
                   _nodesPerCoordinateAxis,_ghostLayerWidth,methodTrace);
      }
    }
  } // Dead code elimination should get rid of this loop if Asserts is not set.
  #endif
}

void exahype::solvers::FiniteVolumesSolver::printFiniteVolumesSolution(
    CellDescription& cellDescription)  const {
  #if DIMENSIONS==2
  double* solution = static_cast<double*>(cellDescription.getSolution());

  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
      std::cout << std::setprecision(3) << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #else
  double* solution = static_cast<double*>(cellDescription.getSolution());

  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
      std::cout << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1 &&
          i(1)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::FiniteVolumesSolver::printFiniteVolumesBoundaryLayer(const double* const luhbnd)  const {
  #if DIMENSIONS==2
  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    for(int i=0; i<_nodesPerCoordinateAxis; ++i) {
      int iScalar = i*_numberOfVariables+unknown;
      std::cout << luhbnd[iScalar] << ",";
      if (i==_nodesPerCoordinateAxis-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #else
  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for(int j=0; j<_nodesPerCoordinateAxis; ++j) {
      for(int i=0; i<_nodesPerCoordinateAxis; ++i) {
        int iScalar = (j*_nodesPerCoordinateAxis+i)*_numberOfVariables+unknown;
        std::cout << luhbnd[iScalar] << ",";
        if (j==_nodesPerCoordinateAxis-1 &&
            i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
      }
    }
    std::cout <<  "}" << std::endl;
  #endif
}

std::string exahype::solvers::FiniteVolumesSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::FiniteVolumesSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:" << _identifier;
  out << ",";
  out << "_type:" << exahype::solvers::Solver::toString(_type);
  out << ",";
  out << "_numberOfVariables:" << _numberOfVariables;
  out << ",";
  out << "_numberOfParameters:" << _numberOfParameters;
  out << ",";
  out << "_nodesPerCoordinateAxis:" << _nodesPerCoordinateAxis;
  out << ",";
  out << "_maximumMeshSize:" << _maximumMeshSize;
  out << ",";
  out << "_timeStepping:" << exahype::solvers::Solver::toString(_timeStepping); // only solver attributes
  out << ",";
  out << "_dataPerPatchFace:" << getDataPerPatchFace();
  out << ",";
  out << "_dataPerPatchBoundary:" << getDataPerPatchBoundary();
  out << ",";
  out << "_dataPerPatch:" << getDataPerPatch();
  out << ",";
  out << "_previousMinTimeStepSize:" << _previousMinTimeStepSize;
  out << ",";
  out << "_minTimeStamp:" << _minTimeStamp;
  out << ",";
  out << "_minTimeStepSize:" << _minTimeStepSize;
  out << ",";
  out << "_nextMinTimeStepSize:" << _admissibleTimeStepSize;
  out <<  ")";
}


void exahype::solvers::FiniteVolumesSolver::putUnknownsIntoByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getExtrapolatedSolutionCompressedIndex()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfExtrapolatedSolution;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> bool  { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getPreviousSolution()),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getExtrapolatedSolution()),
      getDataPerPatchBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
    true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfExtrapolatedSolution);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfExtrapolatedSolution<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
        cellDescription.setPreviousSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion1(
          cellDescription.getPreviousSolutionCompressedIndex()>=0,
          cellDescription.toString()
        );
        lock.free();
        cellDescription.setPreviousSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
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
        cellDescription.setSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion1( cellDescription.getSolutionCompressedIndex()>=0, cellDescription.getSolutionCompressedIndex() );
        lock.free();
        cellDescription.setSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();

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
      cellDescription.setBytesPerDoFInExtrapolatedSolution(compressionOfExtrapolatedSolution);
      if (compressionOfExtrapolatedSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
        cellDescription.setExtrapolatedSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion( cellDescription.getExtrapolatedSolutionCompressedIndex()>=0 );
        lock.free();
        cellDescription.setExtrapolatedSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatchBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedSolutionIndex(), cellDescription.getExtrapolatedSolutionCompressedIndex(), compressionOfExtrapolatedSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolutionIndex(), true );
        cellDescription.setExtrapolatedSolutionIndex(-1);
        cellDescription.setExtrapolatedSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}


void exahype::solvers::FiniteVolumesSolver::pullUnknownsFromByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int dataPerCell     = getDataPerPatch() + getGhostDataPerPatch();
  const int dataPerBoundary = getDataPerPatchBoundary();

  {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell,dataPerCell,DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPerCell,dataPerCell,DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setExtrapolatedSolutionIndex( DataHeap::getInstance().createData( dataPerBoundary,dataPerBoundary,DataHeap::Allocation::UseOnlyRecycledEntries) );
    lock.free();

    if (cellDescription.getPreviousSolutionIndex()==-1) { // allocate new array if recycling has failed
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setExtrapolatedSolutionIndex( DataHeap::getInstance().createData(dataPerBoundary, dataPerBoundary ) );
      lock.free();
    }
    cellDescription.setPreviousSolution    ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionIndex()).data() ) );
    cellDescription.setSolution            ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionIndex()).data() ) );
    cellDescription.setExtrapolatedSolution( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedSolutionIndex()).data() ) );
  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ),
      cellDescription.getPreviousSolutionCompressed()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ),
    cellDescription.getSolutionCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressedIndex() ),
    cellDescription.getExtrapolatedSolutionCompressed()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ), cellDescription.getPreviousSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
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
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
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
      if (cellDescription.getBytesPerDoFInExtrapolatedSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ), cellDescription.getExtrapolatedSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerPatchBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedSolutionIndex(), cellDescription.getExtrapolatedSolutionCompressedIndex(), cellDescription.getBytesPerDoFInExtrapolatedSolution());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolutionCompressedIndex(), true );
        cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
        cellDescription.setExtrapolatedSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	true
  );
}


void exahype::solvers::FiniteVolumesSolver::compress(CellDescription& cellDescription,const bool isSkeletonCell) const {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if ( SpawnCompressionAsBackgroundJob ) {
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);
      peano::datatraversal::TaskSet ( new CompressionJob( *this, cellDescription, isSkeletonCell ));
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(CellDescription::Compressed);
    }
  }
}


void exahype::solvers::FiniteVolumesSolver::uncompress(CellDescription& cellDescription) const {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    madeDecision = cellDescription.getCompressionState() != CellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == CellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( CellDescription::CurrentlyProcessed );
    }
    lock.free();

    peano::datatraversal::TaskSet::finishToProcessBackgroundJobs();
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

  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
}


void exahype::solvers::FiniteVolumesSolver::determineUnknownAverages(
  CellDescription& cellDescription) const {

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionIndex()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionAveragesIndex()), cellDescription.toString() );

  const int dataPerSubcell   = getNumberOfParameters()+getNumberOfVariables();
  const int subcellsPerPatch = (getDataPerPatch()+getGhostDataPerPatch())/ dataPerSubcell;
  const int subcellsPerFace  = (getDataPerPatchFace()) / dataPerSubcell;

  auto& solutionAverages             = getDataHeapEntries( cellDescription.getSolutionAveragesIndex() );
  auto& previousSolutionAverage      = getDataHeapEntries( cellDescription.getPreviousSolutionAveragesIndex() );
  auto& extrapolatedSolutionAverages = getDataHeapEntries( cellDescription.getExtrapolatedSolutionAveragesIndex() );

  // patch data
  kernels::idx2 idx_patchData    (subcellsPerPatch,dataPerSubcell);
  for (int i=0; i<subcellsPerPatch; i++) {
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      solutionAverages[variableNumber]        += getDataHeapEntries(cellDescription.getSolutionIndex())        [idx_patchData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += getDataHeapEntries(cellDescription.getPreviousSolutionIndex())[idx_patchData(i,variableNumber)];
    }
  }
  for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
    solutionAverages[variableNumber]        = solutionAverages[variableNumber]        / (double) subcellsPerPatch;
    previousSolutionAverage[variableNumber] = previousSolutionAverage[variableNumber] / (double) subcellsPerPatch;
  }

  // face data
  kernels::idx2 idx_faceDataAvg(DIMENSIONS_TIMES_TWO,dataPerSubcell);
  kernels::idx3 idx_faceData   (DIMENSIONS_TIMES_TWO,subcellsPerFace,dataPerSubcell);
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<subcellsPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
        extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] +=
            getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex())[idx_faceData(face,i,variableNumber)];
      }
    }
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] =
          extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] / (double) subcellsPerFace;
    }
  }
}


void exahype::solvers::FiniteVolumesSolver::computeHierarchicalTransform(
    CellDescription& cellDescription, double sign) const {
  const int dataPerSubcell   = getNumberOfParameters()+getNumberOfVariables();
  const int subcellsPerPatch = (getDataPerPatch()+getGhostDataPerPatch())/ dataPerSubcell;
  const int subcellsPerFace = (getDataPerPatchFace()) / dataPerSubcell;

  // patch data
  kernels::idx2 idx_patchData    (subcellsPerPatch,dataPerSubcell);
  for (int i=0; i<subcellsPerPatch; i++) {
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      getDataHeapEntries(cellDescription.getSolutionIndex())
          [idx_patchData(i,variableNumber)] += sign * getDataHeapEntries(cellDescription.getSolutionAveragesIndex())[variableNumber];

      getDataHeapEntries(cellDescription.getPreviousSolutionIndex())
          [idx_patchData(i,variableNumber)] += sign * getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex())[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg(DIMENSIONS_TIMES_TWO,dataPerSubcell);
  kernels::idx3 idx_faceData   (DIMENSIONS_TIMES_TWO,subcellsPerFace,dataPerSubcell);
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<subcellsPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ), cellDescription.toString() );
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionAveragesIndex() ), cellDescription.toString() );
        getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex())[idx_faceData(face,i,variableNumber)] +=
              sign * getDataHeapEntries(cellDescription.getExtrapolatedSolutionAveragesIndex())[idx_faceDataAvg(face,variableNumber)];
      }
    }
  }
}

exahype::solvers::FiniteVolumesSolver::CompressionJob::CompressionJob(
  const FiniteVolumesSolver& solver,
  CellDescription&           cellDescription,
  const bool                 isSkeletonJob)
  :
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


bool exahype::solvers::FiniteVolumesSolver::CompressionJob::run(bool runOnMasterThread) {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);
  _cellDescription.setCompressionState(CellDescription::Compressed);

  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_sub(1);
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs.fetch_sub(1);
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
  return false;
}

///////////////////////
// PROFILING
///////////////////////

exahype::solvers::Solver::CellProcessingTimes exahype::solvers::FiniteVolumesSolver::measureCellProcessingTimes(const int numberOfRuns) {
  // Setup
  const int cellDescriptionsIndex = ADERDGSolver::Heap::getInstance().createData(0,1);
  FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionsIndex,0,1);

  Solver::CellInfo cellInfo(cellDescriptionsIndex);
  addNewCellDescription(
      0,cellInfo,CellDescription::Type::Leaf,
      getMaximumAdaptiveMeshLevel(), /* needs to be on the fine grid for the limiter cells */-1,
      getCoarsestMeshSize(),
      _domainOffset);

  CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[0];
  ensureNecessaryMemoryIsAllocated(cellDescription);

  adjustSolutionDuringMeshRefinementBody(cellDescription,true);
  updateTimeStepSize(0,cellInfo);
  cellDescription.setNeighbourMergePerformed(true);

  // MEASUREMENTS
  CellProcessingTimes result;
  const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> boundaryMarkers(0); // >= 0 indicates no remote/domain boundary

  // measure FV cells
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      cellDescription.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      updateBody(cellDescription,cellInfo,true,true,boundaryMarkers,false);

      swapSolutionAndPreviousSolution(cellDescription); // assumed to be very cheap
      rollbackToPreviousTimeStep(cellDescription);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeFVUpdate = time_sec / numberOfRuns;
  }

  // Clean up
  cellDescription.setType(CellDescription::Type::Erased);
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

  DataHeap::getInstance().deleteAllData();
  ADERDGSolver::Heap::getInstance().deleteAllData();
  FiniteVolumesSolver::Heap::getInstance().deleteAllData();

  return result;
}
