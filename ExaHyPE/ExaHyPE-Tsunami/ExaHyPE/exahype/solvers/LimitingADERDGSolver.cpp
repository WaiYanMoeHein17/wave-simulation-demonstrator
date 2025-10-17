/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
  N 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 *
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/

#include <algorithm> // copy_n
#include <iomanip>
#include <chrono>

#include "LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"

#include "kernels/finitevolumes/commons/c/commons.h" // TODO measurements

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

#ifdef USE_ITAC
int exahype::solvers::LimitingADERDGSolver::adjustSolutionHandle            = 0;
int exahype::solvers::LimitingADERDGSolver::fusedTimeStepBodyHandle         = 0;
int exahype::solvers::LimitingADERDGSolver::fusedTimeStepBodyHandleSkeleton = 0;
int exahype::solvers::LimitingADERDGSolver::updateBodyHandle                = 0;
int exahype::solvers::LimitingADERDGSolver::mergeNeighboursHandle           = 0;
#endif

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    exahype::solvers::ADERDGSolver* solver,
    exahype::solvers::FiniteVolumesSolver* limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNumberOfGlobalObservables(), solver->getNodesPerCoordinateAxis(),
        solver->getMaximumMeshSize(),
        solver->getMaximumAdaptiveMeshDepth(),
        solver->getTimeStepping()),
        _solver(std::move(solver)),
        _limiter(std::move(limiter)),
        _DMPMaximumRelaxationParameter(DMPRelaxationParameter),
        _DMPDifferenceScaling(DMPDifferenceScaling)
{
  solver->disableCheckForNaNs();

  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());

}
/** Wire through to limiting ADER-DG solver */
void exahype::solvers::LimitingADERDGSolver::updateMeshUpdateEvent(exahype::solvers::Solver::MeshUpdateEvent meshUpdateEvent) {
  _solver->updateMeshUpdateEvent(meshUpdateEvent);
}
void exahype::solvers::LimitingADERDGSolver::resetMeshUpdateEvent() {
  _solver->resetMeshUpdateEvent();
}
exahype::solvers::Solver::MeshUpdateEvent exahype::solvers::LimitingADERDGSolver::getMeshUpdateEvent() const {
  return _solver->getMeshUpdateEvent();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStamp() const {
  return _solver->getMinTimeStamp();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStepSize() const {
  return _solver->getMinTimeStepSize();
}

double exahype::solvers::LimitingADERDGSolver::getAdmissibleTimeStepSize() const {
  return _solver->getAdmissibleTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::updateAdmissibleTimeStepSize(double value) {
  _solver->updateAdmissibleTimeStepSize(value);
}

void exahype::solvers::LimitingADERDGSolver::resetAdmissibleTimeStepSize() {
  _solver->resetAdmissibleTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::initSolver(
    const double                                timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const double                                boundingBoxSize,
    const double                                boundingBoxMeshSize,
    const std::vector<std::string>&             cmdlineargs,
    const exahype::parser::ParserView&          parserView
) {
  #ifdef Parallel
  // TODO(WORKAROUND): Not sure for what anymore
  const int numberOfObservables = _solver->getDMPObservables();
  _invalidObservables.resize(numberOfObservables);
  std::fill_n(_invalidObservables.data(),_invalidObservables.size(),-1);

  _receivedMax.resize(numberOfObservables);
  _receivedMin.resize(numberOfObservables);
  assertion( numberOfObservables==0 || !_receivedMax.empty());
  assertion( numberOfObservables==0 || !_receivedMin.empty());
  #endif
  
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(
          std::min(boundingBoxMeshSize,_maximumMeshSize),boundingBoxSize);
  _coarsestMeshSize  = coarsestMeshInfo.first; // TODO(Dominic): Wire through as well
  _coarsestMeshLevel = coarsestMeshInfo.second;

  resetMeshUpdateEvent();

  // global observables
  _solver->initSolver(timeStamp, domainOffset, domainSize, boundingBoxSize, boundingBoxMeshSize, cmdlineargs, parserView);
  _limiter->initSolver(timeStamp, domainOffset, domainSize, boundingBoxSize, boundingBoxMeshSize, cmdlineargs, parserView);
}

bool exahype::solvers::LimitingADERDGSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  return _solver->isPerformingPrediction(section);
}

bool exahype::solvers::LimitingADERDGSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  bool isMergingMetadata = false;

  switch (section) {
    case exahype::State::AlgorithmSection::RefinementStatusSpreading:
      isMergingMetadata = getMeshUpdateEvent()==MeshUpdateEvent::IrregularLimiterDomainChange ||
                          getMeshUpdateEvent()==MeshUpdateEvent::RefinementRequested;
      break;
    case exahype::State::AlgorithmSection::MeshRefinement:
      isMergingMetadata = hasRequestedAnyMeshRefinement();
      break;
    default:
      break;
  }

  return isMergingMetadata;
}

void exahype::solvers::LimitingADERDGSolver::kickOffTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch) {
  // copy solver observables to limiter (might be MPI comm. before)
  // TODO add to docu: During mesh refinement, we only require the observables on the main solver
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    std::copy(_solver->_globalObservables.begin(),_solver->_globalObservables.end(),
              _limiter->_globalObservables.begin());
  }

  _solver->kickOffTimeStep(isFirstTimeStepOfBatchOrNoBatch);
  ensureLimiterTimeStepDataIsConsistent();
  _limiter->kickOffTimeStep(isFirstTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::LimitingADERDGSolver::wrapUpTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch,const bool isLastTimeStepOfBatchOrNoBatch) {
  _solver->wrapUpTimeStep(isFirstTimeStepOfBatchOrNoBatch,isLastTimeStepOfBatchOrNoBatch);
  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) { // have consistent data when limiter's wrap up routine is called.
    std::copy(_solver->_globalObservables.begin(),_solver->_globalObservables.end(),
              _limiter->_globalObservables.begin());
  }
  _limiter->wrapUpTimeStep(isFirstTimeStepOfBatchOrNoBatch,isLastTimeStepOfBatchOrNoBatch);
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    SolverPatch& solverPatch,
    Solver::CellInfo& cellInfo) const {
  _solver->synchroniseTimeStepping(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterTimeStepDataIsConsistent() const {
  _limiter->_minTimeStamp            = _solver->_minTimeStamp;
  _limiter->_minTimeStepSize         = _solver->_minTimeStepSize;
  _limiter->_previousMinTimeStamp    = _solver->_previousMinTimeStamp;
  _limiter->_previousMinTimeStepSize = _solver->_previousMinTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSize()  {
  _solver->updateTimeStepSize();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::updateGlobalObservables() {
  _solver->mergeGlobalObservables(_solver->_nextGlobalObservables.data(),_limiter->_nextGlobalObservables.data());
  _solver->updateGlobalObservables();
  std::copy(_solver->_globalObservables.begin(),_solver->_globalObservables.end(),
            _limiter->_globalObservables.begin());
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
  ensureLimiterTimeStepDataIsConsistent();
}

// TODO(Lukas) Is this still needed?
/*
void exahype::solvers::LimitingADERDGSolver::updateNextGlobalObservables(
      const std::vector<double>& globalObservables) {
  Solver::updateNextGlobalObservables(globalObservables);
  _solver->updateNextGlobalObservables(globalObservables);
  // TODO(Lukas) Update for Limiter necessary?
  //_limiter->updateNextGlobalObservables(globalObservables);
}
*/

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatch(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  assertion1(limiterElement!=NotFound,solverPatch.toString());
  LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return cellInfo._FiniteVolumesCellDescriptions[limiterElement];
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatch(
    const SolverPatch& solverPatch,CellInfo& cellInfo,const int limiterElement) const {
  assertion1(limiterElement!=NotFound,solverPatch.toString());
  LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return cellInfo._FiniteVolumesCellDescriptions[limiterElement];
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////

exahype::solvers::Solver::MeshUpdateEvent exahype::solvers::LimitingADERDGSolver::updateRefinementStatusDuringRefinementStatusSpreading(
    SolverPatch& solverPatch) const {
  _solver->updateRefinementStatus(solverPatch);

  if (
       getMeshUpdateEvent()==MeshUpdateEvent::IrregularLimiterDomainChange && // only veto if we expect to a local recomputation, stop if the mesh is not refined
       solverPatch.getLevel() == getMaximumAdaptiveMeshLevel() &&
       ADERDGSolver::isLeaf(solverPatch) &&
       !ADERDGSolver::checkIfStatusFlaggingHasConverged(solverPatch.getFacewiseRefinementStatus(),ADERDGSolver::Erase)
  ) {
    #ifdef MonitorMeshRefinement
    bool converged = ADERDGSolver::checkIfStatusFlaggingHasConverged(solverPatch.getFacewiseRefinementStatus(),ADERDGSolver::Erase);
    logInfo("updateRefinementStatusDuringRefinementStatusSpreading(...)","converged="<<converged);
    if ( !converged ) {
      logInfo("updateRefinementStatusDuringRefinementStatusSpreading(...)","failed for cell="<<solverPatch.toString());
    }
    #endif
    AllSolversAreStable = false;
    if ( tarch::la::min(solverPatch.getFacewiseRefinementStatus())==ADERDGSolver::EmptyStatus ){
      return MeshUpdateEvent::RefinementRequested;
    }
  }

  if ( 
      solverPatch.getType()==SolverPatch::Type::Virtual &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus() > 0
  ) {
    return MeshUpdateEvent::RefinementRequested;
  } else {
    return MeshUpdateEvent::None;
  }
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
 return
     _solver->progressMeshRefinementInEnterCell(
         fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
         coarseGridCell,coarseGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::vetoCoarseningIfRestrictedSolutionIsTroubled(
    const int cellDescriptionsIndex,
    const int fineGridElement) {
  const int DMPObservables = _solver->getDMPObservables();
  if (
      DMPObservables > 0 &&
      fineGridElement!=exahype::solvers::Solver::NotFound
  ) {
    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,fineGridElement);

    if ( solverPatch.getType() == SolverPatch::Type::ParentCoarsens ) {
      auto centre = solverPatch.getOffset() + 0.5 * solverPatch.getSize();
      std::vector<double> DMPObservablesMin(DMPObservables,DMPObservables);
      std::vector<double> DMPObservablesMax(DMPObservables,DMPObservables);

      // previous values
      const double* const previousSolution = static_cast<double*>(solverPatch.getPreviousSolution());
      findCellLocalMinAndMax(previousSolution, DMPObservablesMin.data(), DMPObservablesMax.data());
      bool stable =
          _solver->isPhysicallyAdmissible(previousSolution,DMPObservablesMin.data(),DMPObservablesMax.data(),false,
              centre,solverPatch.getSize(),solverPatch.getPreviousTimeStamp());

      // current values
      const double* const solution = static_cast<double*>(solverPatch.getSolution());
      findCellLocalMinAndMax(previousSolution, DMPObservablesMin.data(), DMPObservablesMax.data());
      stable &=
          _solver->isPhysicallyAdmissible(solution,DMPObservablesMin.data(),DMPObservablesMax.data(),false,
              centre,solverPatch.getSize(),solverPatch.getTimeStamp());

      if ( !stable ) {
        _solver->changeLeafToParent(solverPatch);
        assertion1(solverPatch.getType()==SolverPatch::Type::ParentChecked,solverPatch.toString());
      }
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  const int fineGridElement = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  vetoCoarseningIfRestrictedSolutionIsTroubled(fineGridCell.getCellDescriptionsIndex(),fineGridElement);

  _solver->progressMeshRefinementInLeaveCell(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,fineGridPositionOfCell,solverNumber);
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::LimitingADERDGSolver::eraseOrRefineAdjacentVertices(
    const int cellDescriptionsIndex,
    const int solverNumber,
    const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const bool checkThoroughly,
    bool& checkSuccessful) const {
  return
      _solver->eraseOrRefineAdjacentVertices(
          cellDescriptionsIndex,solverNumber,cellOffset,cellSize,level,checkThoroughly,checkSuccessful);
}

void exahype::solvers::LimitingADERDGSolver::finaliseStateUpdates(
    const int solverNumber,
    CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    if (
        _solver->isLeaf(solverPatch) &&
        solverPatch.getRefinementStatus()<ADERDGSolver::Erase
    ) {
      logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
      std::abort();
    }
    // only done when doing initial refinement
    const bool newLimiterPatchAllocated =
        getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested &&
        ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    if ( newLimiterPatchAllocated ) {
      assertion1(tarch::la::equals(solverPatch.getTimeStamp(),0.0),solverPatch.toString());
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      adjustLimiterSolution(solverPatch,limiterPatch);
    }
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);

    if ( getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested ) {
      solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());
    }

    // must come last because of assertion in patch allocation
    _solver->finaliseStateUpdates(solverNumber,cellInfo);
  }
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////

double exahype::solvers::LimitingADERDGSolver::computeTimeStepSize(SolverPatch& solverPatch,CellInfo& cellInfo) {
  if ( getTimeStepping() == TimeStepping::GlobalFixed ) {
    return solverPatch.getTimeStepSize();
  } else if ( solverPatch.getRefinementStatus() < _solver->_minRefinementStatusForTroubledCell ) {
      return _solver->computeTimeStepSize(solverPatch);
  } else {
    const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
    if ( limiterElement != NotFound ) { // refinement status might have changed just now
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      double* limiterSolution    = static_cast<double*>(limiterPatch.getSolution());
      return _limiter->stableTimeStepSize(limiterSolution, limiterPatch.getSize());
    } else {
      return _solver->computeTimeStepSize(solverPatch);
    }
  }
}

double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(SolverPatch& solverPatch,CellInfo& cellInfo,const bool isFirstTimeStepOfBatch) {
  assertion1(solverPatch.getType()==SolverPatch::Type::Leaf,solverPatch.toString());
  // n-1
  if ( isFirstTimeStepOfBatch ) {
    solverPatch.setPreviousTimeStamp(solverPatch.getTimeStamp());
    solverPatch.setPreviousTimeStepSize(solverPatch.getTimeStepSize());
  }
  // n
  solverPatch.setTimeStamp(solverPatch.getTimeStamp()+solverPatch.getTimeStepSize());

  double admissibleTimeStepSize = computeTimeStepSize(solverPatch,cellInfo);
  solverPatch.setTimeStepSize( admissibleTimeStepSize );
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSize(
    const int solverNumber,
    CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    _solver->updateTimeStepSize(solverNumber,cellInfo);
    ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateGlobalObservables(const int solverNumber,CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    const bool isCell                         = solverPatch.getType()==SolverPatch::Type::Leaf;
    const bool isTroubledCellWithLimiterPatch = isCell &&
        solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell &&
        solverPatch.getPreviousRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2;

    if ( _solver->_numberOfGlobalObservables > 0 && isTroubledCellWithLimiterPatch ) { // TODO(Dominic): Have virtual function; Reassess all that parameter stuff where we need information from user solver; This should just be virtual functions
      const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
      assertion1(limiterElement != NotFound, "Limiter element not found!");
      LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
      const auto cellCentre = solverPatch.getOffset() + 0.5 * solverPatch.getSize();
      const auto& t  = solverPatch.getTimeStamp();
      const auto& dt = solverPatch.getTimeStepSize();
      _limiter->updateGlobalObservables(
          _solver->_nextGlobalObservables.data(),
          static_cast<double*>(limiterPatch.getSolution()),
          cellCentre,solverPatch.getSize(),t,dt); // call must be thread-safe
    } else if ( _solver->_numberOfGlobalObservables > 0 && isCell ) {
      const auto cellCentre = solverPatch.getOffset() + 0.5 * solverPatch.getSize();
      const auto& t  = solverPatch.getTimeStamp();
      const auto& dt = solverPatch.getTimeStepSize();
      _solver->updateGlobalObservables(
          _solver->_nextGlobalObservables.data(),
          static_cast<double*>(solverPatch.getSolution()),
          cellCentre,solverPatch.getSize(),t,dt); // call must be thread-safe
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(SolverPatch& solverPatch,CellInfo& cellInfo) const {
  synchroniseTimeStepping(solverPatch,cellInfo);
  _solver->rollbackToPreviousTimeStep(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinementBody(
    SolverPatch& solverPatch,
    CellInfo& cellInfo,
    const bool isInitialMeshRefinement) {
  #ifdef USE_ITAC
  VT_begin(adjustSolutionHandle);
  #endif

  if (
    solverPatch.getType()==SolverPatch::Type::Leaf ||
    solverPatch.getType()==SolverPatch::Type::LeafProlongates
  ) {
    _solver->ensureNecessaryMemoryIsAllocated(solverPatch);
    if (solverPatch.getType()==SolverPatch::Type::LeafProlongates) {
      _solver->prolongateVolumeData(solverPatch,isInitialMeshRefinement);
    }

    _solver->adjustSolution(solverPatch);
    determineSolverMinAndMax(solverPatch,false);
    if ( !evaluatePhysicalAdmissibilityCriterion(solverPatch,solverPatch.getTimeStamp()) ) {
       solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
    } else {
      _solver->markForRefinement(solverPatch);
    }
    solverPatch.setType(SolverPatch::Type::LeafChecked);
  }

  solverPatch.setHasCompletedLastStep(true);

  #ifdef USE_ITAC
  VT_end(adjustSolutionHandle);
  #endif
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if ( limiterElement!=Solver::NotFound ) {
    LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
    copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  }
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) {
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, LimiterPatch& limiterPatch) {
  limiterPatch.setPreviousTimeStamp(solverPatch.getPreviousTimeStamp());
  limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousTimeStepSize());
  limiterPatch.setTimeStamp(solverPatch.getTimeStamp());
  limiterPatch.setTimeStepSize(solverPatch.getTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::fusedTimeStepOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int element        = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(solverPatch,cellInfo);
    solverPatch.setHasCompletedLastStep(false);

    // Write the previous limiter status back onto the patch for all cell description types
    solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

    if ( ADERDGSolver::isLeaf(solverPatch) ) {
      const bool isAMRSkeletonCell     = solverPatch.getAugmentationStatus() > ADERDGSolver::MinimumAugmentationStatusForVirtualRefining;
      const bool isAtRemoteBoundary    = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
      const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
      const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

      if (
          (SpawnUpdateAsBackgroundJob || (SpawnPredictionAsBackgroundJob && !isLastTimeStepOfBatch)) &&
          !mustBeDoneImmediately
      ) {
        const auto predictionTimeStepData = _solver->getPredictionTimeStepData(solverPatch,true/*duringFusedTimeStep*/);
        peano::datatraversal::TaskSet( new FusedTimeStepJob(
            *this,solverPatch,cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch,isLastTimeStepOfBatch,boundaryMarkers,isSkeletonCell) );
      } else {
        const auto predictionTimeStepData = _solver->getPredictionTimeStepData(solverPatch,true/*duringFusedTimeStep*/);
        fusedTimeStepBody(
            solverPatch, cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch,isLastTimeStepOfBatch,
            boundaryMarkers,isSkeletonCell,
            mustBeDoneImmediately);
      }
    }
    else { // other cell types
      if (
          solverPatch.getType()==SolverPatch::Type::Virtual &&
          solverPatch.getCommunicationStatus()>=ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication
      ) {
        _solver->restrictToTopMostParent(solverPatch,isFirstTimeStepOfBatch/*addToCoarseGridUpdate*/);
      }
      ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
      updateMeshUpdateEvent( _solver->updateRefinementStatusAfterSolutionUpdate(solverPatch) );
      solverPatch.setHasCompletedLastStep(true);
    }
  }
}

/**
 * @return Comptues a merged limiter status as a maximum of
 * the current cell value and the neighbours' values decremented by 1.
 *
 * @param[in] solverPatch A solver patch.
 */
int exahype::solvers::LimitingADERDGSolver::computeMergedRefinementStatus(const SolverPatch& solverPatch) {
  int newRefinementStatus = solverPatch.getRefinementStatus();
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    newRefinementStatus = std::max(newRefinementStatus,solverPatch.getFacewiseRefinementStatus(faceIndex)-1);
  }
  return newRefinementStatus;
}

void exahype::solvers::LimitingADERDGSolver::fusedTimeStepBody(
    SolverPatch&                                       solverPatch,
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

  updateSolution(solverPatch,cellInfo,isFirstTimeStepOfBatch,boundaryMarkers,
                 isFirstTimeStepOfBatch/*addSurfaceIntegralContributionToUpdate*/);
  const bool isTroubled = checkIfCellIsTroubledAndDetermineMinAndMax(solverPatch,cellInfo);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(solverPatch,cellInfo,isFirstTimeStepOfBatch);
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(solverPatch,cellInfo,isTroubled);

  reduce(solverPatch,cellInfo,result);

  if (
      solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell &&
      SpawnPredictionAsBackgroundJob &&
      !mustBeDoneImmediately         &&
      isLastTimeStepOfBatch  // may only spawned in last iteration
  ) {
    const int element = cellInfo.indexOfADERDGCellDescription(solverPatch.getSolverNumber());
    peano::datatraversal::TaskSet(
        new ADERDGSolver::PredictionJob(
            *_solver.get(),solverPatch/*the reductions are delegated to _solver anyway*/,
            cellInfo._cellDescriptionsIndex,element,
            predictionTimeStamp,
            predictionTimeStepSize,
            false/*is uncompressed*/,isSkeletonCell,isLastTimeStepOfBatch/*addVolumeIntegralResultToUpdate*/));
  }
  else if ( solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell ){
    _solver->predictionAndVolumeIntegralBody(
        solverPatch,
        predictionTimeStamp,
        predictionTimeStepSize,
        false/*is uncompressed*/,isSkeletonCell,isLastTimeStepOfBatch/*addVolumeIntegralResultToUpdate*/);
  } else {
    solverPatch.setHasCompletedLastStep(true);
  }

  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_end(fusedTimeStepBodyHandleSkeleton);
  } else {
    VT_end(fusedTimeStepBodyHandle);
  }
  #endif
}

void exahype::solvers::LimitingADERDGSolver::predictionAndVolumeIntegral(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(solverPatch,cellInfo);

    const bool isAMRSkeletonCell = ADERDGSolver::belongsToAMRSkeleton(solverPatch);
    const bool isSkeletonCell    = isAMRSkeletonCell || isAtRemoteBoundary;
    waitUntilCompletedLastStep(solverPatch,isSkeletonCell,false);
    if ( _solver->isLeaf(solverPatch) ) {
      if ( solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell ) { // only compute predictor for cells which need to communicate with ADER-DG neighbours
        const auto predictionTimeStepData = _solver->getPredictionTimeStepData(solverPatch,false); // this is either the fused scheme or a predictor recomputation
        _solver->predictionAndVolumeIntegral(
            solverNumber,cellInfo,
            std::get<0>(predictionTimeStepData),
            std::get<1>(predictionTimeStepData),
            true,isAtRemoteBoundary,
            FuseAllADERDGPhases || areRollbacksPossible()/*addVolumeIntegralResultToUpdate*/);
      }
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::resetGlobalObservables(
    double* const globalObservables) {
  // refers call to the wrapped solvers
  _solver->resetGlobalObservables();
  _limiter->resetGlobalObservables();
}

void exahype::solvers::LimitingADERDGSolver::updateGlobalObservables(
    double* const                               globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const double t,
    const double dt) {
  logError("resetGlobalObservables(...)","routine never be called!");
  std::abort();
}

void exahype::solvers::LimitingADERDGSolver::mergeGlobalObservables(
    double* const       globalObservables,
    const double* const otherObservables) {
  logError("mergeGlobalObservables(...)","routine should never be called!");
  std::abort();
}

void exahype::solvers::LimitingADERDGSolver::wrapUpGlobalObservables(
    double* const globalObservables) {
  logError("wrapUpGlobalObservables(...)","routine should never be called!");
  std::abort();
}

void exahype::solvers::LimitingADERDGSolver::reduce(
    const SolverPatch&  solverPatch,
    CellInfo&           cellInfo,
    const UpdateResult& result) {
  // all reductions are written to the solver's global fields
  // the limiter's fields are synchronised in kickOffTimeStep
  updateAdmissibleTimeStepSize(result._timeStepSize);
  updateMeshUpdateEvent(result._meshUpdateEvent);
  updateGlobalObservables(solverPatch.getSolverNumber(),cellInfo);
}

void exahype::solvers::LimitingADERDGSolver::updateBody(
    SolverPatch&                                       solverPatch,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers){
  #ifdef USE_ITAC
  VT_begin(updateBodyHandle);
  #endif

  if (CompressionAccuracy>0.0) { uncompress(solverPatch,cellInfo); }

  // the actual computations
  updateSolution(solverPatch,cellInfo,true,boundaryMarkers,areRollbacksPossible()/*effect: add surface integral result to solution*/);
  const bool isTroubled = checkIfCellIsTroubledAndDetermineMinAndMax(solverPatch,cellInfo);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(solverPatch,cellInfo,true); // uses DG solution to compute time step size; might be result of FV->DG projection
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(solverPatch,cellInfo,isTroubled);

  reduce(solverPatch,cellInfo,result);

  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
  if (CompressionAccuracy>0.0) { compress(solverPatch,cellInfo,isAtRemoteBoundary); }

  solverPatch.setHasCompletedLastStep(true); // required as prediction checks the flag too. Field should be renamed "setHasCompletedLastOperation(...)".

  #ifdef USE_ITAC
  VT_end(updateBodyHandle);
  #endif
}

void exahype::solvers::LimitingADERDGSolver::updateOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    synchroniseTimeStepping(solverPatch,cellInfo);
    solverPatch.setHasCompletedLastStep(false);

    const bool isAtRemoteBoundary    = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
    if ( _solver->isLeaf(solverPatch) && SpawnUpdateAsBackgroundJob ) {
      peano::datatraversal::TaskSet( new UpdateJob(*this, solverPatch,cellInfo,boundaryMarkers ) );
    }
    else if ( _solver->isLeaf(solverPatch) ) {
      updateBody(solverPatch,cellInfo,boundaryMarkers);
    }
    else { // other cell types
      if (
          solverPatch.getType()==SolverPatch::Type::Virtual &&
          solverPatch.getCommunicationStatus()>=ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication
      ) {
        _solver->restrictToTopMostParent(solverPatch,areRollbacksPossible()/*effect: add surface integral result to solution*/);
      }
      ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
      updateMeshUpdateEvent(_solver->updateRefinementStatusAfterSolutionUpdate(solverPatch) );
      solverPatch.setHasCompletedLastStep(true);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::uncompress(
    SolverPatch& solverPatch,CellInfo& cellInfo) const {
  _solver->uncompress(solverPatch);
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
    _limiter->uncompress(limiterPatch);
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    compress(solverPatch,cellInfo,isAtRemoteBoundary);
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    SolverPatch& solverPatch,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  if ( _solver->isLeaf(solverPatch) ) {
    _solver->compress(solverPatch,isAtRemoteBoundary);
    const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
    if ( limiterElement!=exahype::solvers::Solver::NotFound ) {
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
      _limiter->compress(limiterPatch,isAtRemoteBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinement(
    const int solverNumber,
    CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    synchroniseTimeStepping(solverPatch,cellInfo);

    const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
    if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
      solverPatch.setHasCompletedLastStep(false);
      peano::datatraversal::TaskSet(
        new AdjustSolutionDuringMeshRefinementJob(*this,solverPatch,cellInfo,isInitialMeshRefinement) );
    } else {
      adjustSolutionDuringMeshRefinementBody(solverPatch,cellInfo,isInitialMeshRefinement);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::updateSolution(
    SolverPatch&                                       solverPatch,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStep,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         addSurfaceIntegralResultToSolution) {
  assertion(solverPatch.getRefinementStatus()>=ADERDGSolver::Erase);
  const int& limiterStatus = solverPatch.getRefinementStatus();
  const bool isTroubledCellOrDirectNeighbour =
      _solver->isLeaf(solverPatch) &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      limiterStatus >= _solver->_minRefinementStatusForTroubledCell-1;
  const bool isSecondDegreeNeighbourOfTroubled =
      _solver->isLeaf(solverPatch) &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      limiterStatus == _solver->_minRefinementStatusForTroubledCell-2;

  if ( isTroubledCellOrDirectNeighbour ) { // limiter update
    assertion1(solverPatch.getRefinementStatus()>0,solverPatch.toString());
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);

    _limiter->updateSolution(limiterPatch,cellInfo._cellDescriptionsIndex,boundaryMarkers,isFirstTimeStep);

    // always perform the projection in neighbours of troubled cells as these couple the dg and fv domain.
    // only if both are deactivated, no rollbacks are performed and the troubled cells can skip the projection.
    const bool isNeighourOfTroubledOrRollbacksRequired =
        !(OnlyStaticLimiting && OnlyInitialMeshRefinement) ||
        limiterStatus == _solver->_minRefinementStatusForTroubledCell-1;
    if ( isNeighourOfTroubledOrRollbacksRequired ) {
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    }
  } else {
    _solver->correction(solverPatch,boundaryMarkers,isFirstTimeStep,addSurfaceIntegralResultToSolution);

    if ( isSecondDegreeNeighbourOfTroubled ) {
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch); // will do rollback if troubled in separation layers
    } else {
      ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
    }
  }
}

bool
exahype::solvers::LimitingADERDGSolver::checkIfCellIsTroubledAndDetermineMinAndMax(
    SolverPatch& solverPatch,
    CellInfo&    cellInfo) {
  if ( OnlyStaticLimiting ) {
    return solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell;
  }

  bool isTroubled =
      !evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) ||
      !evaluatePhysicalAdmissibilityCriterion(solverPatch,
                                              solverPatch.getTimeStamp()+solverPatch.getTimeStepSize()); // after min and max was found

  if ( // above call computes DG min and max on-the-fly. We use the FV min and max if the cell is troubled
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell
  ) {
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    determineLimiterMinAndMax(solverPatch,limiterPatch);
  }

  return isTroubled;
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::updateRefinementStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,
    CellInfo&    cellInfo,
    const bool   isTroubled) {
  assertion1(_solver->isLeaf(solverPatch),solverPatch.toString());

  if ( OnlyInitialMeshRefinement && OnlyStaticLimiting ) {
    return MeshUpdateEvent::None;
  }

  // pre-update mesh update events
  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
  if ( isTroubled &&
       solverPatch.getLevel() == getMaximumAdaptiveMeshLevel() &&
       solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell
  ) {
    meshUpdateEvent = MeshUpdateEvent::IrregularLimiterDomainChange;
  } else if (
      isTroubled &&
      solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
  ) {
    meshUpdateEvent = MeshUpdateEvent::RefinementRequested;
  }

  // update refinement status
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

  _solver->updateRefinementStatus(solverPatch);
  if ( isTroubled ) {
    solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
  }
  else if ( // cool down
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getPreviousRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell
  ) {
    solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell-2); // allows increasing the status if we are "suddenly" next to troubled cell
    solverPatch.setRefinementStatus(computeMergedRefinementStatus(solverPatch));
  }

  // post-update mesh update events (pure ADER-DG functionality
  if (
      solverPatch.getLevel()<getMaximumAdaptiveMeshLevel() ||
      solverPatch.getRefinementStatus() < _solver->_minRefinementStatusForTroubledCell-2
  ) {
    meshUpdateEvent = Solver::mergeMeshUpdateEvents(
        meshUpdateEvent,
        _solver->updateRefinementStatusAfterSolutionUpdate(solverPatch));
  }

  return meshUpdateEvent;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = static_cast<double*>(solverPatch.getSolution());

  const int numberOfObservables = _solver->getDMPObservables();
  if ( numberOfObservables>0 ) {
    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = discreteMaximumPrincipleAndMinAndMaxSearch(solution, observablesMin,observablesMax);

    // 1. Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    // todo later on we might evaluate the DMP also during the mesh refinement iterations.
    // Then, we might need to pass the time stamp as well

    return
        dmpIsSatisfied ||
        _solver->vetoDiscreteMaximumPrincipleDecision(
            solution,
            observablesMin,observablesMax,
            solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell,
            solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
            solverPatch.getTimeStamp()+solverPatch.getTimeStepSize());
  } else {
    return true;
  }
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(
    SolverPatch& solverPatch,const double timeStamp) {
  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    observablesMax = static_cast<double*>(solverPatch.getSolutionMax());
  }

  const double* const solution = static_cast<double*>(solverPatch.getSolution());

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      timeStamp);
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int solverNumber,Solver::CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    if (
        _solver->isLeaf(solverPatch) &&
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
    ) {
      assertion1( solverPatch.getRefinementStatus()>=ADERDGSolver::Erase,solverPatch.getRefinementStatus() );

      if ( solverPatch.getRefinementStatus()<ADERDGSolver::Erase ) {
        logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
        std::abort();
      } else if (solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell) {
        determineSolverMinAndMax(solverPatch,true);
      } else {
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
    } else if ( _solver->isLeaf(solverPatch) ) {
      determineSolverMinAndMax(solverPatch,true);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(
    SolverPatch& solverPatch, const bool validate) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionIndex()),
            solverPatch.toString());
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMinIndex()),
            solverPatch.toString());

    const double* const solution = static_cast<double*>(solverPatch.getSolution());

    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

    // Write the result to the face with index "0"
    findCellLocalMinAndMax(solution, observablesMin, observablesMax);

    // Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    #ifdef Asserts
    if ( validate ) {
      for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
        assertion2(*(observablesMin+i)<std::numeric_limits<double>::infinity(),i,solverPatch.toString());
        assertion2(*(observablesMax+i)>-std::numeric_limits<double>::infinity(),i,solverPatch.toString());
      } // Dead code elimination will get rid of this loop
    }
    #endif
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* limiterSolution = static_cast<double*>(limiterPatch.getSolution());

    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

    // Write the result to the face with index "0"
    findCellLocalLimiterMinAndMax(limiterSolution, observablesMin, observablesMax);

    // Copy the result on the other faces as well
    const int numberOfObservables = _solver->getDMPObservables();
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion(*(observablesMin+i)<std::numeric_limits<double>::infinity());
      assertion(*(observablesMax+i)>-std::numeric_limits<double>::infinity());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatch(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  assertion(limiterElement!=Solver::NotFound);
  LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
  limiterPatch.setType(LimiterPatch::Type::Erased);
  _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch); // has its own lock
  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  cellInfo._FiniteVolumesCellDescriptions.erase(cellInfo._FiniteVolumesCellDescriptions.begin()+limiterElement);
  lock.free();
}

void exahype::solvers::LimitingADERDGSolver::ensureNoLimiterPatchIsAllocatedOnHelperCell(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement!=exahype::solvers::Solver::NotFound &&
      !_solver->isLeaf(solverPatch)
  ) {
    deallocateLimiterPatch(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement!=Solver::NotFound      &&
      solverPatch.getRefinementStatus()         < _solver->_minRefinementStatusForTroubledCell-2 &&
      solverPatch.getPreviousRefinementStatus() < _solver->_minRefinementStatusForTroubledCell-2
  ) {
    deallocateLimiterPatch(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustLimiterSolution(
    SolverPatch& solverPatch,
    LimiterPatch& limiterPatch) const {
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  _limiter->adjustSolution(limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  #if defined(Asserts)
  const int previousLimiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  #endif
  assertion(previousLimiterElement==Solver::NotFound);

  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      solverPatch.getSolverNumber(),
      cellInfo,
      LimiterPatch::Type::Leaf,
      solverPatch.getLevel(),
      solverPatch.getParentIndex(),
      solverPatch.getSize(),
      solverPatch.getOffset());

  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolutionIndex()),solverPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionIndex()),solverPatch.toString());

  LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

  assertion1(DataHeap::getInstance().isValidIndex(limiterPatch.getPreviousSolutionIndex()),limiterPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(limiterPatch.getSolutionIndex()),limiterPatch.toString());
  assertion1(limiterPatch.getSolution()!=nullptr,limiterPatch.toString());
  assertion1(limiterPatch.getPreviousSolution()!=nullptr,limiterPatch.toString());
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const SolverPatch& solverPatch,
        CellInfo&          cellInfo,
        const int          limiterStatus) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement        ==Solver::NotFound              &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      _solver->isLeaf(solverPatch)                          &&
      limiterStatus         >= _solver->_minRefinementStatusForTroubledCell-2
  ) {
    assertion1(solverPatch.getPreviousRefinementStatus()<_solver->_minRefinementStatusForTroubledCell-2,
               solverPatch.toString());
    allocateLimiterPatch(solverPatch,cellInfo);
    return true;
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = static_cast<double*>(solverPatch.getSolution());
  double*       limiterSolution = static_cast<double*>(limiterPatch.getSolution());

  projectOnFVLimiterSpace(solverSolution, limiterSolution);
}

void exahype::solvers::LimitingADERDGSolver::rollbackSolutionGlobally(const int solverNumber,CellInfo& cellInfo) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    rollbackToPreviousTimeStep(solverPatch,cellInfo);

    // 1. Ensure limiter patch is allocated (based on previous limiter status
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getPreviousRefinementStatus());

    // 2. Roll solution back to previous time step
    if ( ADERDGSolver::isLeaf(solverPatch) ) {
      if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
        assertion1(_solver->isLeaf(solverPatch),solverPatch.toString());

        _solver->swapSolutionAndPreviousSolution(solverPatch);   // roll back solver

        if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-2 ) {
          LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
          _limiter->swapSolutionAndPreviousSolution(limiterPatch); // roll back limiter (must exist!)
        } else {
          const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
          if ( limiterElement!=Solver::NotFound ) {
            LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
            copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
            projectDGSolutionOnFVSpace(solverPatch,limiterPatch); // project DG solution on patch
          }
        }
      }
      else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
        _solver->swapSolutionAndPreviousSolution(solverPatch);
      }
    }

    // 3. Reset the current to the previous refinement status.
    solverPatch.setRefinementStatus(solverPatch.getPreviousRefinementStatus());
  }
}

void exahype::solvers::LimitingADERDGSolver::rollbackSolutionLocally(
    const int  solverNumber,
    CellInfo&  cellInfo,
    const bool fusedTimeStepping) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    // 1. Ensure limiter patch is allocated (based on current limiter status)
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());

    // 2. Now roll back to the last valid solution
    if ( isInvolvedInLocalRecomputation(solverPatch) ) {  // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
      assertion1(_solver->isLeaf(solverPatch),solverPatch.toString());
      if ( !_solver->isLeaf(solverPatch) ) {
        logError("rollbackSolutionLocally(..)","type is not cell but is=" << solverPatch.toString());
        std::abort();
      }

      if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-1 ) {
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      }
      else { // We need to project limiter data for the previous stage
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
    }

    // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = static_cast<double*>(limiterPatch.getSolution());
  double*       solverSolution  = static_cast<double*>(solverPatch.getSolution());

  projectOnDGSpace(limiterSolution, solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::localRecomputation(
    SolverPatch&                                       solverPatch,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  assertion1(isInvolvedInLocalRecomputation(solverPatch),solverPatch.toString());
  if ( solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-1 ) {
    // these guys are recomputing with the limiter
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);

    _limiter->rollbackToPreviousTimeStep(limiterPatch);
    _limiter->updateSolution(limiterPatch,cellInfo._cellDescriptionsIndex,boundaryMarkers,true);
    copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch); // restore pre-rollback limiter patch time stamp and step size
    projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
  }
  else {
    // these guys are just swapping and projecting
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    _solver->swapSolutionAndPreviousSolution(solverPatch);
    _limiter->swapSolutionAndPreviousSolution(limiterPatch);
    projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
  }
}

double exahype::solvers::LimitingADERDGSolver::localRecomputationBody(
    SolverPatch&                                       solverPatch,
    Solver::CellInfo&                                  cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  // 1. Perform the local recomputation in affected cells
  if ( isInvolvedInLocalRecomputation(solverPatch) ) {
    localRecomputation(solverPatch,cellInfo,boundaryMarkers);
  }

  // 2. Re-compute a new time step size in ALL Cells
  double admissibleTimeStepSize = std::numeric_limits<double>::infinity();
  if ( _solver->isLeaf(solverPatch) ) {
    admissibleTimeStepSize = computeTimeStepSize(solverPatch,cellInfo);
    solverPatch.setTimeStepSize(admissibleTimeStepSize);
    ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
  }

  // 3. Recompute the predictor in certain cells if fused time stepping is used.
  const bool isNeigbourOfTroubledOrWasPreviouslyTroubled =
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      _solver->isLeaf(solverPatch) &&
      solverPatch.getRefinementStatus() < _solver->_minRefinementStatusForTroubledCell &&
      (solverPatch.getRefinementStatus()        == _solver->_minRefinementStatusForTroubledCell-1 ||
      solverPatch.getPreviousRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell);
  if (
      FuseAllADERDGPhases &&
      isNeigbourOfTroubledOrWasPreviouslyTroubled
  ) {
    const auto predictionTimeStepData = _solver->getPredictionTimeStepData(solverPatch,true/*duringFusedTimeStep*/);
    const bool isAtRemoteBoundary     = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
    _solver->predictionAndVolumeIntegral(
        solverPatch.getSolverNumber(),cellInfo,
        std::get<0>(predictionTimeStepData),
        std::get<1>(predictionTimeStepData),
        false/*already uncompressed*/,isAtRemoteBoundary,true/*addVolumeIntegralResultToUpdate*/);
  }

  // 4. Determine min and max in recomputed cellss
  if ( solverPatch.getRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-1 ) {
    determineMinAndMax(solverPatch.getSolverNumber(),cellInfo);
  }

  // 5. Set the previous limiter status to the current one
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());
  return admissibleTimeStepSize;
}


void exahype::solvers::LimitingADERDGSolver::localRecomputation(
    const int                                          solverNumber,
    Solver::CellInfo&                                  cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int solverElement  = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    // 1. Perform the local recomputation
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    if ( isInvolvedInLocalRecomputation(solverPatch) ) {
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    }
    if ( SpawnUpdateAsBackgroundJob ) {
      peano::datatraversal::TaskSet( new LocalRecomputationJob(
          *this, solverPatch,cellInfo,boundaryMarkers ) );
    } else {
      double admissibleTimeStepSize = localRecomputationBody(
          solverPatch,cellInfo,boundaryMarkers);
      updateAdmissibleTimeStepSize(admissibleTimeStepSize);
    }
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursDataDuringLocalRecomputation(
    const int                                  solverNumber,
    Solver::CellInfo&                          cellInfo1,
    Solver::CellInfo&                          cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>&  pos1,
    const tarch::la::Vector<DIMENSIONS, int>&  pos2) {
  const int solverElement1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int solverElement2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);

  if ( solverElement1!=Solver::NotFound && solverElement2!=Solver::NotFound) {
    SolverPatch& solverPatch1 = cellInfo1._ADERDGCellDescriptions[solverElement1];
    SolverPatch& solverPatch2 = cellInfo2._ADERDGCellDescriptions[solverElement2];

    if (_solver->isLeaf(solverPatch1) &&
        _solver->isLeaf(solverPatch2) &&
        solverPatch1.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2 &&
        solverPatch2.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2) {
      assertion2(solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());
      assertion2(solverPatch2.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());

      waitUntilCompletedLastStep<SolverPatch>(solverPatch1,false,false); // must come before any other operation
      waitUntilCompletedLastStep<SolverPatch>(solverPatch2,false,false);

      const int minStatus = std::min(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      const int maxStatus = std::max(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      if ( minStatus<=_solver->_minRefinementStatusForTroubledCell-2 &&
           maxStatus>=_solver->_minRefinementStatusForTroubledCell ) {
        logError("mergeNeighboursDataDuringLocalRecomputation(...)","Neighbours cannot communicate during local recomputation." <<
            std::endl << "cell1=" << solverPatch1.toString() <<
            std::endl << ".cell2=" << solverPatch2.toString());
        std::terminate();
      }

      _limiter->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighboursData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo1,
    Solver::CellInfo&                         cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  #ifdef USE_ITAC
  VT_begin(mergeNeighboursHandle);
  #endif

  const int solverElement1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int solverElement2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);

  if ( solverElement1!=Solver::NotFound && solverElement2!=Solver::NotFound) {
    SolverPatch& solverPatch1 = cellInfo1._ADERDGCellDescriptions[solverElement1];
    SolverPatch& solverPatch2 = cellInfo2._ADERDGCellDescriptions[solverElement2];

    waitUntilCompletedLastStep<SolverPatch>(solverPatch1,false,false); // must come before any other operation
    waitUntilCompletedLastStep<SolverPatch>(solverPatch2,false,false);

    //
    // 1. Solve the Riemann problems/copy boundary layers
    //
    // We only limit on the finest mesh level
    if ( solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion2(solverPatch2.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());

      if (solverPatch1.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell &&
          solverPatch2.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell) {
        // assumes that face fluxes are not directly added to cell update or the update is cleared if FV update is performe
        _solver->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
      }

      if (solverPatch1.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2 &&
          solverPatch2.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2) {
        _limiter->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
      }

      const int minStatus = std::min(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      const int maxStatus = std::max(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      if ( minStatus<=_solver->_minRefinementStatusForTroubledCell-2 &&
           maxStatus>=_solver->_minRefinementStatusForTroubledCell ) {
        logError("mergeNeighboursBasedOnLimiterStatus(...)","Neighbours cannot communicate. " <<
            std::endl << "cell1=" << solverPatch1.toString() <<
            std::endl << ".cell2=" << solverPatch2.toString());
        std::terminate();
      }
    // On the other levels, we work with the ADER-DG solver only
    } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      _solver->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
    }

    // 2. Merge the min and max of both cell description's solver's
    // solution value. This is not done during recomputations.
    if ( _solver->getDMPObservables() > 0 ) {
      Solver::InterfaceInfo face(pos1,pos2);
      mergeSolutionMinMaxOnFace(solverPatch1,solverPatch2,face);
    }
  }

  #ifdef USE_ITAC
  VT_end(mergeNeighboursHandle);
  #endif
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  Solver::InterfaceInfo& face) const {
  if ( ADERDGSolver::communicateWithNeighbour(solverPatch1,face._faceIndex1) ) {
    assertion1( ADERDGSolver::communicateWithNeighbour(solverPatch2,face._faceIndex2),solverPatch2.toString() );
    assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
    const int numberOfObservables = _solver->getDMPObservables();
    double* min1 = static_cast<double*>(solverPatch1.getSolutionMin()) + face._faceIndex1 * numberOfObservables;
    double* min2 = static_cast<double*>(solverPatch2.getSolutionMin()) + face._faceIndex2 * numberOfObservables;
    double* max1 = static_cast<double*>(solverPatch1.getSolutionMax()) + face._faceIndex1 * numberOfObservables;
    double* max2 = static_cast<double*>(solverPatch2.getSolutionMax()) + face._faceIndex2 * numberOfObservables;

    for (int i=0; i<numberOfObservables; i++) {
      const double min = std::min(
          *(min1+i),
          *(min2+i)
      );
      const double max = std::max(
          *(max1+i),
          *(max2+i)
      );

      *(min1+i) = min;
      *(min2+i) = min;

      *(max1+i) = max;
      *(max2+i) = max;
    }
  }
}

#ifdef Parallel
///////////////////////////////////
// NEIGHBOUR - Mesh refinement
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  _solver->appendNeighbourCommunicationMetadata(
      metadata,src,dest,
      cellDescriptionsIndex,solverNumber);
}

/////////////////////////////////////
// NEIGHBOUR - Local Recomputation
/////////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourDuringLocalRecomputation(
        const int                                    toRank,
        const int                                    solverNumber,
        CellInfo&                                    cellInfo,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const tarch::la::Vector<DIMENSIONS, double>& barycentre,
        const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(src,dest);

  if (
      level==getMaximumAdaptiveMeshLevel() &&
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    logDebug("sendDataToNeighbourDuringLocalRecomputation(...)", "send data for solver " << _identifier << " to rank="<<toRank<<",x="<<barycentre<<",level="<<level);

    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    waitUntilCompletedLastStep<SolverPatch>(solverPatch,true,true);

    if ( solverPatch.getRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-2 ) {
      _limiter->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
    } else {
      _limiter->sendEmptyDataToNeighbour(toRank,barycentre,level);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataDuringLocalRecomputation(
    const int                                    fromRank,
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(dest,src); // ! order of arguments

  if (
      level==getMaximumAdaptiveMeshLevel() &&
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    logDebug("mergeWithNeighbourDataDuringLocalRecomputation(...)", "receive data for solver " << _identifier << " from rank="<<fromRank<<",x="<<barycentre<<",level="<<level);

    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    if ( solverPatch.getRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-1 ) {
      _limiter->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
    } else {
      _limiter->dropNeighbourData(fromRank,barycentre,level);
    }
  }
}

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                    toRank,
    const SolverPatch&                           solverPatch,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  BoundaryFaceInfo face(src,dest);
  const int numberOfObservables = _solver->getDMPObservables();
  if (
      numberOfObservables>0 &&
      ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex)
  ) {
    assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMinIndex()));
    assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMaxIndex()));
    const double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin()) + (face._faceIndex * numberOfObservables);
    const double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax()) + (face._faceIndex * numberOfObservables);

    DataHeap::getInstance().sendData(
        observablesMin, numberOfObservables, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().sendData(
        observablesMax, numberOfObservables, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                    toRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const int                                    level) {
  BoundaryFaceInfo face(src,dest);
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    // wait; must come before any other operation
    waitUntilCompletedLastStep<SolverPatch>(solverPatch,true,true);

    // send order:   minAndMax,solver,limiter
    // receive order: limiter,solver,minAndMax
    // 1. send min and max
    sendMinAndMaxToNeighbour(toRank,solverPatch,src,dest,barycentre,level);
    // 2. solver sends
    _solver->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
    // 3. limiter sends (receive order must be inverted)
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank="<<toRank<<",x="<<barycentre<<",level="<<level);

      if (
          solverPatch.getRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-2 &&
          cellInfo.indexOfFiniteVolumesCellDescription(solverNumber) != Solver::NotFound // might not be allocated yet
      ) {
        _limiter->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
      } else {
        _limiter->sendEmptyDataToNeighbour(toRank,barycentre,level);
      }
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  solverPatch,
  Solver::BoundaryFaceInfo& face,
  const double* const min,
  const double* const max) const {
  assertion1(ADERDGSolver::ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex) ,solverPatch.toString());

  double* solutionMin = static_cast<double*>(solverPatch.getSolutionMin());
  double* solutionMax = static_cast<double*>(solverPatch.getSolutionMax());

  const int numberOfObservables = _solver->getDMPObservables();
  for (int i=0; i<numberOfObservables; i++) {
    solutionMin[i+face._faceIndex*numberOfObservables]  = std::min( solutionMin[i+face._faceIndex*numberOfObservables], min[i] );
    solutionMax[i+face._faceIndex*numberOfObservables]  = std::max( solutionMax[i+face._faceIndex*numberOfObservables], max[i] );
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    SolverPatch&                                 solverPatch,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  BoundaryFaceInfo face(dest,src);
  if (
      numberOfObservables>0 &&
      ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex)
  ) {
    // Inverted send-receive order: TODO(Dominic): Add to docu
    // Send order:    min,max
    // Receive order; max,min
    DataHeap::getInstance().receiveData(
        const_cast<double*>(_receivedMax.data()), numberOfObservables, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(
        const_cast<double*>(_receivedMin.data()), numberOfObservables, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

    mergeSolutionMinMaxOnFace(solverPatch,face,_receivedMin.data(),_receivedMax.data());
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  logDebug("mergeWithNeighbourData(...)", "receive for solver " << _identifier <<" from rank="<<fromRank<<"x="<<x<<",level="<<level);
  BoundaryFaceInfo face(dest,src);
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    assertion1(solverPatch.getRefinementStatus()>=ADERDGSolver::Pending,solverPatch.toString());

    // send order:   minAndMax,solver,limiter
    // receive order limiter,solver,minAndMax

    // 1. limiter receives
    if ( level == getMaximumAdaptiveMeshLevel() ) {
      if (
          solverPatch.getRefinementStatus()                        >= _solver->_minRefinementStatusForTroubledCell-2 &&
          solverPatch.getFacewiseRefinementStatus(face._faceIndex) >= _solver->_minRefinementStatusForTroubledCell-2
      ) {
        assertion1(cellInfo.indexOfFiniteVolumesCellDescription(solverNumber)!=Solver::NotFound,solverPatch.toString());
        _limiter->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
      } else {
        _limiter->dropNeighbourData(fromRank,x,level);
      }
    }

    // 2. solver receives
    if (
        level < getMaximumAdaptiveMeshLevel()
        ||
        (solverPatch.getRefinementStatus() < _solver->_minRefinementStatusForTroubledCell &&
        solverPatch.getFacewiseRefinementStatus(face._faceIndex) < _solver->_minRefinementStatusForTroubledCell)
    ) {
      _solver->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
    } else {
      _solver->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
    }

    // 3. min and max
    mergeWithNeighbourMinAndMax(fromRank,cellInfo._ADERDGCellDescriptions[solverElement],src,dest,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(dest,src);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    // send order:   minAndMax,solver,limiter
    // receive order limiter,solver,minAndMax
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      _limiter->dropNeighbourData(fromRank,x,level);
    }
    _solver->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);

    const int numberOfObservables = _solver->getDMPObservables();
    if ( numberOfObservables>0 ) {
      for(int receives=0; receives<2; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

/////////////////////////////////////
// MASTER<=>WORKER
/////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  _solver->progressMeshRefinementInPrepareSendToWorker(
      workerRank, fineGridCell, fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell, coarseGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerIfProlongating(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToWorkerIfProlongating(
      workerRank,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::receiveDataFromMasterIfProlongating(
    const int masterRank,
    const int receivedCellDescriptionsIndex,
    const int receivedElement,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  _solver->receiveDataFromMasterIfProlongating(
      masterRank,receivedCellDescriptionsIndex,receivedElement,x,level);
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
 return _solver->progressMeshRefinementInMergeWithWorker(
      localCellDescriptionsIndex,
      receivedCellDescriptionsIndex,receivedElement);
}

void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  _solver->progressMeshRefinementInPrepareSendToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int  level) {
  vetoCoarseningIfRestrictedSolutionIsTroubled(localCellDescriptionsIndex,localElement);

  _solver->progressMeshRefinementInMergeWithMaster(
      worker,localCellDescriptionsIndex,localElement,coarseGridCellDescriptionsIndex,x,level);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,messageType,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
        toRank,cellDescriptionsIndex,limiterElement,messageType,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,messageType,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,messageType,x,level);
  }
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  _solver->sendDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  _solver->sendDataToWorker(workerRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithMasterData(masterRank,x,level);
}
#endif

std::string exahype::solvers::LimitingADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::LimitingADERDGSolver::toString (std::ostream& out) const {
  out << getIdentifier() << "{_ADERDG: ";
  out << _solver->toString() << "}\n";
  out << getIdentifier() << "{_FV: ";
  out << _limiter->toString() << "}";
}

exahype::solvers::Solver::CellProcessingTimes exahype::solvers::LimitingADERDGSolver::measureCellProcessingTimes(const int numberOfRuns) {
  // Setup
  const int cellDescriptionsIndex = ADERDGSolver::Heap::getInstance().createData(0,1);
  FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionsIndex,0,1);

  Solver::CellInfo cellInfo(cellDescriptionsIndex);
  _solver->addNewCellDescription(
      0,cellInfo,SolverPatch::Type::Leaf,
      getMaximumAdaptiveMeshLevel(), /* needs to be on the fine grid for the limiter cells */-1,
      getCoarsestMeshSize(),
      _domainOffset);

  SolverPatch& solverPatch   = cellInfo._ADERDGCellDescriptions[0];
  _solver->ensureNecessaryMemoryIsAllocated(solverPatch);

  adjustSolutionDuringMeshRefinementBody(solverPatch,cellInfo,true);
  updateTimeStepSize(0,cellInfo);
  const double dt = _solver->_admissibleTimeStepSize;
  _solver->_minTimeStepSize       = dt;
  _solver->_estimatedTimeStepSize = dt;

  // ADER-DG specific setup ( all Riemanns have been performed, cell is surrounded by other Cell type cells )
  solverPatch.setNeighbourMergePerformed(true);
  solverPatch.setAugmentationStatus(0);
  solverPatch.setFacewiseAugmentationStatus(0);
  solverPatch.setCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);
  solverPatch.setFacewiseCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);

  // MEASUREMENTS
  CellProcessingTimes result;
  const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> boundaryMarkers(0); // >= 0 indicates no remote/domain boundary

  // measure ADERDG STP
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    int numberOfPicardIterations = std::numeric_limits<int>::max();
    for (int it=0; it<numberOfRuns; it++) {
      numberOfPicardIterations = _solver->predictionAndVolumeIntegralBody(
          solverPatch,0,dt,false,true,true);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._minTimePredictor = time_sec / numberOfRuns / std::abs(numberOfPicardIterations);
    if ( _solver->isLinear() ) {
      result._maxTimePredictor = result._minTimePredictor;
    } else {
      result._maxTimePredictor = result._minTimePredictor * getNodesPerCoordinateAxis(); // * (order+1)
    }
  }

  // measure ADER-DG cells
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_refineOrKeepOnFineGrid);
      solverPatch.setFacewiseRefinementStatus(_solver->_refineOrKeepOnFineGrid); // assumed  to be very cheap

      solverPatch.setTimeStamp(0);
      solverPatch.setTimeStepSize(dt);
      solverPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      updateBody(solverPatch,cellInfo,boundaryMarkers);

      _solver->swapSolutionAndPreviousSolution(solverPatch); // assumed  to be very cheap
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure ADER-DG -> FV cells
  solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell-2);
  solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell-2);
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    solverPatch.setTimeStamp(0);
    solverPatch.setTimeStepSize(dt);
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> centre(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(centre,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-centre);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_refineOrKeepOnFineGrid);
      solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell-2); // assumed  to be very cheap

      solverPatch.setTimeStamp(0);
      solverPatch.setTimeStepSize(dt);
      solverPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      limiterPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      updateBody(solverPatch,cellInfo,boundaryMarkers);

      _solver->swapSolutionAndPreviousSolution(solverPatch); // assumed  to be very cheap
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDG2FVUpdate = time_sec / numberOfRuns;
  }

  // measure FV -> ADERDG cells
  solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell-1);
  solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    solverPatch.setTimeStamp(0);
    solverPatch.setTimeStepSize(dt);
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> centre(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(centre,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-centre);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell-1);
      solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell); // assumed  to be very cheap

      solverPatch.setTimeStamp(0);
      solverPatch.setTimeStepSize(dt);
      solverPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      limiterPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      updateBody(solverPatch,cellInfo,boundaryMarkers);

      _solver->swapSolutionAndPreviousSolution(solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch); // assumed  to be very cheap
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeFV2ADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure troubled / FV cells
  solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
  solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    solverPatch.setTimeStamp(0);
    solverPatch.setTimeStepSize(dt);
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> centre(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(centre,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-centre);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
      solverPatch.setFacewiseRefinementStatus(_solver->_minRefinementStatusForTroubledCell);

      solverPatch.setTimeStamp(0);
      solverPatch.setTimeStepSize(dt);
      solverPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      limiterPatch.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      updateBody(solverPatch,cellInfo,boundaryMarkers);

      _solver->swapSolutionAndPreviousSolution(solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeFVUpdate = time_sec / numberOfRuns;
  }

  // measure Riemann solve
  {
    const tarch::la::Vector<DIMENSIONS,int> pos1(0);
    tarch::la::Vector<DIMENSIONS,int> pos2(0); pos2[0]=1;
    Solver::InterfaceInfo face(pos1,pos2);
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      _solver->solveRiemannProblemAtInterface(solverPatch,solverPatch,face);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGRiemann = time_sec / numberOfRuns;
  }

  // Clean up
  solverPatch.setRefinementStatus(0);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
  solverPatch.setType(SolverPatch::Type::Erased);
  _solver->ensureNoUnnecessaryMemoryIsAllocated(solverPatch);

  DataHeap::getInstance().deleteAllData();
  ADERDGSolver::Heap::getInstance().deleteAllData();
  FiniteVolumesSolver::Heap::getInstance().deleteAllData();

  return result;
}
