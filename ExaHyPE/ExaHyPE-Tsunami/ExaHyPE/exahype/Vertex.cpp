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
 
#include "exahype/Vertex.h"

#include "tarch/la/Scalar.h"

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/dForRange.h"

#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"

#include "exahype/State.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"

#include <iterator>

tarch::logging::Log exahype::Vertex::_log( "exahype::Vertex");

exahype::Vertex::Vertex() : Base() {
  _vertexData.setCellDescriptionsIndex(mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex);
  _vertexData.setADERDGCellDescriptions(static_cast<void*>(nullptr));
  _vertexData.setFiniteVolumesCellDescriptions(static_cast<void*>(nullptr));
}

exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
  // do nothing
}

exahype::Vertex::Vertex(const Base::PersistentVertex& argument)
    : Base(argument) {
  // do nothing
}

bool exahype::Vertex::SpawnNeighbourMergeAsThread = false;

peano::MappingSpecification exahype::Vertex::getNeighbourMergeSpecification(const int level) {
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel ) {
    return peano::MappingSpecification(
           peano::MappingSpecification::WholeTree,
           peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
}

bool exahype::Vertex::equalUpToRelativeTolerance(
    const tarch::la::Vector<DIMENSIONS,double>& lhs,
    const tarch::la::Vector<DIMENSIONS,double>& rhs) {
  const double tolerance =
      tarch::la::NUMERICAL_ZERO_DIFFERENCE *
      std::max(
          1.0, std::max( tarch::la::maxAbs(lhs), tarch::la::maxAbs(rhs) )
  );
  return tarch::la::equals( lhs, rhs, tolerance );
}

tarch::la::Vector<TWO_POWER_D, int> exahype::Vertex::getCellDescriptionsIndex() const {
  return _vertexData.getCellDescriptionsIndex();
}

int exahype::Vertex::getCellDescriptionsIndex(const int adjacencyIndex) const {
  return _vertexData.getCellDescriptionsIndex(adjacencyIndex);
}

exahype::solvers::Solver::CellInfo exahype::Vertex::createCellInfo(int index) const {
  if ( _vertexData.getCellDescriptionsIndex(index)<0 ) {
    logError("createCellInfo()","No heap data allocated for cell. Cell descriptions heap index="<<_vertexData.getCellDescriptionsIndex(index))
    std::abort();
  }
  assertion1( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_vertexData.getCellDescriptionsIndex(index)),_vertexData.getCellDescriptionsIndex());
  assertion2( (index >=0) && (index < TWO_POWER_D), index, TWO_POWER_D );

  return solvers::Solver::CellInfo(
      _vertexData.getCellDescriptionsIndex(index),_vertexData.getADERDGCellDescriptions(index),_vertexData.getFiniteVolumesCellDescriptions(index));
}

tarch::la::Vector<DIMENSIONS,double> exahype::Vertex::computeFaceBarycentre(
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const tarch::la::Vector<DIMENSIONS,double>& h,
    const int                                   normalDirection,
    const tarch::la::Vector<DIMENSIONS,int>&    cellPosition) {
  tarch::la::Vector<DIMENSIONS,double> barycentre;
  for (int d=0; d<DIMENSIONS; d++) {
    barycentre(d) = x(d) + 0.5*h(d)*(2*cellPosition(d)-1);
  }
  barycentre(normalDirection) = x(normalDirection);
  return barycentre;
}

exahype::solvers::Solver::RefinementControl exahype::Vertex::evaluateRefinementCriterion(
    const tarch::la::Vector<DIMENSIONS, double>& vertexOffset,
    const int                                    level,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const bool                                   checkThoroughly,
    bool&                                        checkSuccessful) const {
  bool canErase   = true;
  bool mustRefine = false;
  dfor2(pos)
    tarch::la::Vector<DIMENSIONS,double> cellOffset(vertexOffset);
    for (int d = 0; d < DIMENSIONS; ++d) {
      cellOffset[d] += ( pos[d]-1 ) * cellSize[d];
    }

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      const int cellDescriptionsIndex = _vertexData.getCellDescriptionsIndex(posScalar);
      bool geometryInfoCorrect = false;
      exahype::solvers::Solver::RefinementControl control =
          solver->eraseOrRefineAdjacentVertices(
              cellDescriptionsIndex,solverNumber,cellOffset,cellSize,level,checkThoroughly,geometryInfoCorrect);
      mustRefine      |= (control==exahype::solvers::Solver::RefinementControl::Refine);
      canErase        &= (control==exahype::solvers::Solver::RefinementControl::Erase);
      checkSuccessful &= geometryInfoCorrect;
    }
  enddforx

  if (mustRefine) {
    return exahype::solvers::Solver::RefinementControl::Refine;
  } else if (canErase) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  } else {
    return exahype::solvers::Solver::RefinementControl::Keep;
  }
}

void exahype::Vertex::mergeOnlyNeighboursMetadataLoopBodyHelper(
    solvers::Solver::CellInfo&                   cellInfo1,
    const solvers::Solver::CellInfo&             cellInfo2,
    const tarch::la::Vector<DIMENSIONS,int>&     pos1,
    const tarch::la::Vector<DIMENSIONS,int>&     pos2,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const exahype::State::AlgorithmSection&      section) {
  for ( auto& patch1 : cellInfo1._ADERDGCellDescriptions ) {
    auto* solver1 = solvers::RegisteredSolvers[patch1.getSolverNumber()];
    if ( solver1->isMergingMetadata(section) ) {
      for ( auto& patch2 : cellInfo2._ADERDGCellDescriptions ) {
        solvers::ADERDGSolver::mergeWithNeighbourMetadata(patch1.getSolverNumber(),cellInfo1,
            patch2.getAugmentationStatus(),patch2.getCommunicationStatus(),patch2.getRefinementStatus(), // patch2 values
            pos1,pos2,barycentre);
      }
    }
  }
}

void exahype::Vertex::mergeOnlyNeighboursMetadataLoopBody(
    const int pos1Scalar,
    const int pos2Scalar,
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const exahype::State::AlgorithmSection& section,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const bool                                   checkThoroughly) {
  tarch::la::Vector<DIMENSIONS,int> pos1 = delineariseIndex2(pos1Scalar);
  tarch::la::Vector<DIMENSIONS,int> pos2 = delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  const bool validIndex1 = cellDescriptionsIndex1 >= 0;
  const bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    solvers::Solver::CellInfo cellInfo1(cellDescriptionsIndex1);
    solvers::Solver::CellInfo cellInfo2(cellDescriptionsIndex2);

    solvers::Solver::InterfaceInfo face(pos1,pos2);
    const auto barycentre = computeFaceBarycentre(x,h,face._direction,pos1);
    mergeOnlyNeighboursMetadataLoopBodyHelper(cellInfo1,cellInfo2,pos1,pos2,barycentre,section);
    mergeOnlyNeighboursMetadataLoopBodyHelper(cellInfo2,cellInfo1,pos2,pos1,barycentre,section);
  } else if (
      validIndex1 != validIndex2 &&
      cellDescriptionsIndex1 != mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex &&
      cellDescriptionsIndex2 != mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex
  ) {
    const int cellDescriptionsIndex      = validIndex1 ? cellDescriptionsIndex1 : cellDescriptionsIndex2;
    const int otherCellDescriptionsIndex = validIndex1 ? cellDescriptionsIndex2 : cellDescriptionsIndex1;
    tarch::la::Vector<DIMENSIONS,int> pos      = validIndex1 ? pos1 : pos2;
    tarch::la::Vector<DIMENSIONS,int> posEmpty = validIndex1 ? pos2 : pos1;

    solvers::Solver::CellInfo cellInfo(cellDescriptionsIndex);
    for ( auto& patch : cellInfo._ADERDGCellDescriptions) {
      auto* solver = solvers::RegisteredSolvers[patch.getSolverNumber()];
      if ( solver->isMergingMetadata(section) ) {
        solvers::Solver::InterfaceInfo face(pos1,pos2);
        const auto barycentre = computeFaceBarycentre(x,h,face._direction,pos1);
        solvers::ADERDGSolver::mergeWithEmptyNeighbourDuringMeshRefinement(
            patch.getSolverNumber(),cellInfo,pos,posEmpty,
            otherCellDescriptionsIndex==mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex,
            barycentre);
      }
    }
  }
}

void exahype::Vertex::mergeOnlyNeighboursMetadata(
    const exahype::State::AlgorithmSection& section,
    const tarch::la::Vector<DIMENSIONS,     double>& x,
    const tarch::la::Vector<DIMENSIONS,     double>& h,
    const bool                              checkThoroughly) const {
  assertion(!isHangingNode());
  assertion(isInside() || isBoundary());
  mergeOnlyNeighboursMetadataLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),section,x,h,checkThoroughly);
  #if DIMENSIONS==3
  mergeOnlyNeighboursMetadataLoopBody(0,4,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(4),section,x,h,checkThoroughly);
  #endif
}

void exahype::Vertex::validateNeighbourhood(
    const int                                cellDescriptionsIndex1,
    const int                                cellDescriptionsIndex2,
    const exahype::Vertex&                   vertex,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2) {
  tarch::la::Vector<DIMENSIONS,int> posCell  = pos1;
  tarch::la::Vector<DIMENSIONS,int> posEmpty = pos2;
  if ( cellDescriptionsIndex2 >= 0 ) {
    posCell  = pos2;
    posEmpty = pos1;
  }
  const int posCellScalar = peano::utils::dLinearised(posCell,2);
  solvers::Solver::CellInfo cellInfo = vertex.createCellInfo(posCellScalar);
  solvers::Solver::BoundaryFaceInfo face(posCell,posEmpty);

  // ADER-DG
  for (auto& p : cellInfo._ADERDGCellDescriptions) {
    if (
        (solvers::ADERDGSolver::isLeaf(p) ||
         solvers::ADERDGSolver::isParent(p))
        &&
        cellDescriptionsIndex2==mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex
    ) {
      logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex<<" next to empty cell: cell="<<p.toString());
      std::terminate();
    }
  }
  // FV
  for (auto& p : cellInfo._FiniteVolumesCellDescriptions) {
    if (
        p.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Leaf &&
        cellDescriptionsIndex2==mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex
    ) {
      logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex<<" next to empty cell: cell="<<p.toString());
      std::terminate();
    }
  }
}

void exahype::Vertex::mergeNeighboursDataAndMetadata(
    solvers::Solver::CellInfo& cellInfo1,
    solvers::Solver::CellInfo& cellInfo2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  if ( !cellInfo1.empty() && !cellInfo2.empty() ) {
    solvers::Solver::InterfaceInfo face(pos1,pos2);
    const auto barycentre = computeFaceBarycentre(x,h,face._direction,pos1);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          static_cast<solvers::ADERDGSolver*>(solver)->
              mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,barycentre,false);
          static_cast<solvers::ADERDGSolver*>(solver)->
              mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
          break;
        case solvers::Solver::Type::LimitingADERDG:
          static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
              mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,barycentre,false);
          static_cast<solvers::LimitingADERDGSolver*>(solver)->
              mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
          break;
        case solvers::Solver::Type::FiniteVolumes:
          static_cast<solvers::FiniteVolumesSolver*>(solver)->
          mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
          break;
        default:
          assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          logError("mergeNeighboursDataAndMetadata(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          std::abort();
          break;
      }
    }
  }
}

tarch::la::Vector<DIMENSIONS,int> exahype::Vertex::delineariseIndex2(int index) {
  #if DIMENSIONS==2
  switch( index ) {                                //x,y
    case 0: return tarch::la::Vector<DIMENSIONS,int>(0,0);
    case 1: return tarch::la::Vector<DIMENSIONS,int>(1,0);
    case 2: return tarch::la::Vector<DIMENSIONS,int>(0,1);
    case 3: return tarch::la::Vector<DIMENSIONS,int>(1,1);
    default:
      logError("delineariseIndex2(index)","index must be in range [0,3]!");
      return tarch::la::Vector<DIMENSIONS,int>(-1,-1);
  }
  #elif DIMENSIONS==3
  switch( index ) {                                //x,y,z
    case 0: return tarch::la::Vector<DIMENSIONS,int>(0,0,0);
    case 1: return tarch::la::Vector<DIMENSIONS,int>(1,0,0);
    case 2: return tarch::la::Vector<DIMENSIONS,int>(0,1,0);
    case 3: return tarch::la::Vector<DIMENSIONS,int>(1,1,0);
    case 4: return tarch::la::Vector<DIMENSIONS,int>(0,0,1);
    case 5: return tarch::la::Vector<DIMENSIONS,int>(1,0,1);
    case 6: return tarch::la::Vector<DIMENSIONS,int>(0,1,1);
    case 7: return tarch::la::Vector<DIMENSIONS,int>(1,1,1);
    default:
      logError("delineariseIndex2(index)","index must be in range [0,7]!");
      return tarch::la::Vector<DIMENSIONS,int>(-1,-1,-1);
  }
  #endif
}

void exahype::Vertex::mergeNeighboursLoopBody(
    const int                                    pos1Scalar,
    const int                                    pos2Scalar,
    const exahype::Vertex&                       vertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  #if defined(Asserts) || defined (ValidateNeighbourHoodDuringNeighbourMerge)
  constexpr bool validate = true;
  #else
  constexpr bool validate = false;
  #endif

  const int cellDescriptionsIndex1 = VertexOperations::readCellDescriptionsIndex(vertex,pos1Scalar);
  const int cellDescriptionsIndex2 = VertexOperations::readCellDescriptionsIndex(vertex,pos2Scalar);
  assertion2(cellDescriptionsIndex1 < 1 || cellDescriptionsIndex1 != cellDescriptionsIndex2,cellDescriptionsIndex1,cellDescriptionsIndex2);

  const tarch::la::Vector<DIMENSIONS,int> pos1 = delineariseIndex2(pos1Scalar);
  const tarch::la::Vector<DIMENSIONS,int> pos2 = delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  bool validIndex1 = cellDescriptionsIndex1 >= 0;
  bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    solvers::Solver::CellInfo cellInfo1 = vertex.createCellInfo(pos1Scalar);
    solvers::Solver::CellInfo cellInfo2 = vertex.createCellInfo(pos2Scalar);
    mergeNeighboursDataAndMetadata(cellInfo1,cellInfo2,pos1,pos2,x,h);
  }
  else if  (
      validate
      &&
      validIndex1 != validIndex2
      &&
      cellDescriptionsIndex1!=mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex &&
      cellDescriptionsIndex2!=mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex
      &&
      cellDescriptionsIndex1!=mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex &&
      cellDescriptionsIndex2!=mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex
  ) {
    validateNeighbourhood(cellDescriptionsIndex1,cellDescriptionsIndex2,vertex,pos1,pos2);
  }
}

void exahype::Vertex::mergeNeighbours(
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  if ( tarch::la::allSmallerEquals(h,exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers()) ) {
    if ( SpawnNeighbourMergeAsThread ) {
      #if DIMENSIONS==2
      peano::datatraversal::TaskSet runParallelTasks(
      [&]() -> bool { mergeNeighboursLoopBody(0,1,*this,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,2,*this,x,h); return false; },
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      true);
      #elif DIMENSIONS==3
      peano::datatraversal::TaskSet runParallelTasks(
      [&]() -> bool { mergeNeighboursLoopBody(0,1,*this,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,2,*this,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,4,*this,x,h); return false; },
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      true);
      #endif
    } else {
     mergeNeighboursLoopBody(0,1,*this,x,h);
     mergeNeighboursLoopBody(0,2,*this,x,h);
     #if DIMENSIONS==3
     mergeNeighboursLoopBody(0,4,*this,x,h);
     #endif
    }
  }
}

// PARALLEL

#if Parallel
bool exahype::Vertex::hasToCommunicate( const int level) const {
  return (level >= exahype::solvers::Solver::getCoarsestMeshLevelOfAllSolvers());
}

bool exahype::Vertex::hasToSendMetadata(
  const int                                 toRank,
  const int                                 srcScalar,
  const int                                 destScalar,
  const tarch::la::Vector<TWO_POWER_D,int>& adjacentRanks) {
  return adjacentRanks[destScalar] == toRank
         &&
         adjacentRanks[destScalar] != tarch::parallel::Node::getGlobalMasterRank() &&
         adjacentRanks[srcScalar]  != tarch::parallel::Node::getGlobalMasterRank()
         &&
         (   // Send also when a fork/join was triggered for the current rank
             adjacentRanks[srcScalar]   == tarch::parallel::Node::getInstance().getRank() ||
             State::isForkTriggeredForRank(adjacentRanks[srcScalar]) ||
             State::isJoinTriggeredForRank(adjacentRanks[srcScalar])
         );
}

bool exahype::Vertex::compareGeometryInformationOfCellDescriptionsAndVertex(
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  if ( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) ) {
    solvers::Solver::CellInfo cellInfo(srcCellDescriptionsIndex);

    if ( !cellInfo.empty() ) {
      solvers::Solver::BoundaryFaceInfo face(src,dest);

      bool geometryInformationMatches = true;

      for ( auto& patch : cellInfo._ADERDGCellDescriptions ) {
        tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
            exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
        geometryInformationMatches &= equalUpToRelativeTolerance(barycentreFromPatch,barycentre);
      }
      for ( auto& patch : cellInfo._FiniteVolumesCellDescriptions ) {
        tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
            exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
        geometryInformationMatches &= equalUpToRelativeTolerance(barycentreFromPatch,barycentre);
      }
      return geometryInformationMatches;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

void exahype::Vertex::sendOnlyMetadataToNeighbourLoopBody(
    const int                                    toRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    srcCellDescriptionsIndex,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const bool                                   checkThoroughly) {
  const auto src  = delineariseIndex2(srcScalar);
  const auto dest = delineariseIndex2(destScalar);

  if ( hasToSendMetadata(toRank,srcScalar,destScalar,adjacentRanks) ) {
    solvers::Solver::BoundaryFaceInfo face(dest,src);
    const auto barycentre = computeFaceBarycentre(x,h,face._direction,dest);

    logDebug("sendOnlyMetadataToNeighbourLoopBody(...)","to rank="<<toRank<<",barycentre="<<barycentre.toString()<<",level="<<level<<",adjacentRanks="<<adjacentRanks);

    if ( !checkThoroughly || compareGeometryInformationOfCellDescriptionsAndVertex(src,dest,srcCellDescriptionsIndex,barycentre,h) ) {
      exahype::sendNeighbourCommunicationMetadata(toRank,srcCellDescriptionsIndex,src,dest,barycentre,level);
    } else {
      exahype::sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(toRank,barycentre,level);
    }
  }
}

void exahype::Vertex::sendOnlyMetadataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const bool                                   checkThoroughly) const {
  if ( hasToCommunicate(level) ) {
    const tarch::la::Vector<DIMENSIONS,int> lowerLeft(0);
    sendOnlyMetadataToNeighbourLoopBody(toRank,0,1,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,checkThoroughly);
    sendOnlyMetadataToNeighbourLoopBody(toRank,1,0,_vertexData.getCellDescriptionsIndex(1),getAdjacentRanks(),x,h,level,checkThoroughly);
    sendOnlyMetadataToNeighbourLoopBody(toRank,0,2,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,checkThoroughly);
    sendOnlyMetadataToNeighbourLoopBody(toRank,2,0,_vertexData.getCellDescriptionsIndex(2),getAdjacentRanks(),x,h,level,checkThoroughly);
    #if DIMENSIONS==3
    sendOnlyMetadataToNeighbourLoopBody(toRank,0,4,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,checkThoroughly);
    sendOnlyMetadataToNeighbourLoopBody(toRank,4,0,_vertexData.getCellDescriptionsIndex(4),getAdjacentRanks(),x,h,level,checkThoroughly);
    #endif
  }
}


bool exahype::Vertex::hasToReceiveMetadata(
    const int                                  fromRank,
    const int                                  srcScalar,
    const int                                  destScalar,
    const tarch::la::Vector<TWO_POWER_D, int>& adjacentRanks) {
  return
      adjacentRanks[srcScalar]    == fromRank
      &&
      adjacentRanks[srcScalar]    != tarch::parallel::Node::getGlobalMasterRank() &&
      adjacentRanks[destScalar]   != tarch::parallel::Node::getGlobalMasterRank()
      &&
      ( // Receive also when a fork/join is performed for the neighbour rank, i.e. such
        // an event was triggered in the iteration before
        adjacentRanks[destScalar] == tarch::parallel::Node::getInstance().getRank() ||
        State::isForkingRank(adjacentRanks[destScalar]) ||
        State::isJoiningRank(adjacentRanks[destScalar])
      );
}


void exahype::Vertex::mergeOnlyWithNeighbourMetadataLoopBody(
    const int                                    fromRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const exahype::State::AlgorithmSection&      section,
    const bool                                   checkThoroughly) {
  const auto src  = delineariseIndex2(srcScalar);
  const auto dest = delineariseIndex2(destScalar);

  if ( hasToReceiveMetadata(fromRank,srcScalar,destScalar,adjacentRanks) ) {
    solvers::Solver::BoundaryFaceInfo face(dest,src);
    const auto barycentre = computeFaceBarycentre(x,h,face._direction,dest);

    logDebug("mergeOnlyWithNeighbourMetadataLoopBody(...)","from rank="<<fromRank<<",barycentre="<<barycentre.toString()<<",level="<<level<<",adjacentRanks="<<adjacentRanks);

    exahype::MetadataHeap::HeapEntries receivedMetadata;
    receivedMetadata.clear();
    exahype::receiveNeighbourCommunicationMetadata(receivedMetadata,fromRank,barycentre,level);
    assertionEquals(receivedMetadata.size(),exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

    bool validIndex = destCellDescriptionIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(destCellDescriptionIndex);

      for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
        auto* solver = solvers::RegisteredSolvers[solverNumber];
        if ( solver->isMergingMetadata(section) ) {
          const int begin = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
          const int end   = begin+exahype::NeighbourCommunicationMetadataPerSolver;
          exahype::MetadataHeap::HeapEntries metadataPortion(receivedMetadata.begin()+begin,receivedMetadata.begin()+end);

          switch ( solver->getType() ) {
            case solvers::Solver::Type::ADERDG: // break  removed by purpose
            case solvers::Solver::Type::LimitingADERDG: {
              solvers::ADERDGSolver::mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,dest,src,barycentre);
            } break;
            case solvers::Solver::Type::FiniteVolumes:
              // do nothing
              break;
            default:
              assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              std::abort();
              break;
          }
        }
      }
    }
  }
}

void exahype::Vertex::mergeOnlyWithNeighbourMetadata(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const exahype::State::AlgorithmSection&      section,
    const bool                                   checkThoroughly) const {
  if ( hasToCommunicate(level) ) {
    #if DIMENSIONS==3
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,4,0,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,0,4,_vertexData.getCellDescriptionsIndex(4),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    #endif
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,2,0,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,0,2,_vertexData.getCellDescriptionsIndex(2),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,1,0,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    mergeOnlyWithNeighbourMetadataLoopBody(fromRank,0,1,_vertexData.getCellDescriptionsIndex(1),getAdjacentRanks(),x,h,level,section,checkThoroughly);
  }
}

void exahype::Vertex::dropNeighbourMetadata(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    const tarch::la::Vector<DIMENSIONS,int> lowerLeft(0); // dest and src is swapped & order is swapped
    #if DIMENSIONS==3
    if ( hasToReceiveMetadata(fromRank,4,0,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,2,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
    if ( hasToReceiveMetadata(fromRank,0,4,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,2,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
    #endif
    if ( hasToReceiveMetadata(fromRank,2,0,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,1,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
    if ( hasToReceiveMetadata(fromRank,0,2,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,1,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
    if ( hasToReceiveMetadata(fromRank,1,0,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,0,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
    if ( hasToReceiveMetadata(fromRank,0,1,getAdjacentRanks()) ) { MetadataHeap::getInstance().receiveData(fromRank,computeFaceBarycentre(x,h,0,lowerLeft),level,peano::heap::MessageType::NeighbourCommunication); }
  }
}

void exahype::Vertex::sendToNeighbourLoopBody(
  const int                                    toRank,
  const int                                    srcScalar,
  const int                                    destScalar,
  const int                                    srcCellDescriptionsIndex,
  const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
  const bool                                   isLastIterationOfBatchOrNoBatch,
  const tarch::la::Vector<DIMENSIONS, double>& x,
  const tarch::la::Vector<DIMENSIONS, double>& h,
  const int                                    level) {
  const auto src  = delineariseIndex2(srcScalar);
  const auto dest = delineariseIndex2(destScalar);

  solvers::Solver::BoundaryFaceInfo face(src,dest);
  const auto barycentre = computeFaceBarycentre(x,h,face._direction,src);

  if ( hasToSendMetadata(toRank,srcScalar,destScalar,adjacentRanks) ) {
    logDebug("sendToNeighbourLoopBody(...)","to rank="<<toRank<<",barycentre="<<barycentre.toString()<<",src="<<src.toString()<<",dest="<<dest.toString());

    const bool validIndex = srcCellDescriptionsIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex));

    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(srcCellDescriptionsIndex);
      for ( unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber ) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        switch ( solver->getType() ) {
          case solvers::Solver::Type::ADERDG:
            static_cast<solvers::ADERDGSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
            break;
          case solvers::Solver::Type::LimitingADERDG:
            static_cast<solvers::LimitingADERDGSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
            break;
          case solvers::Solver::Type::FiniteVolumes:
            static_cast<solvers::FiniteVolumesSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,barycentre,level);
            break;
          default:
            assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
            logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
            std::abort();
            break;
        }
      }
    }

    // metadata is sent and received as block
    const bool sendMetadata = !solvers::Solver::DisableMetadataExchangeDuringTimeSteps;
    if ( sendMetadata && validIndex ){
      sendNeighbourCommunicationMetadata(
          toRank,srcCellDescriptionsIndex,src,dest,barycentre,level);
    } else if ( sendMetadata ) {
      sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(toRank,barycentre,level);
    }
  }
}

void exahype::Vertex::sendToNeighbour(
    int                                          toRank,
    bool                                         isLastIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    sendToNeighbourLoopBody(toRank,0,1,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    sendToNeighbourLoopBody(toRank,1,0,_vertexData.getCellDescriptionsIndex(1),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    sendToNeighbourLoopBody(toRank,0,2,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    sendToNeighbourLoopBody(toRank,2,0,_vertexData.getCellDescriptionsIndex(2),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    #if DIMENSIONS==3
    sendToNeighbourLoopBody(toRank,0,4,_vertexData.getCellDescriptionsIndex(0),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    sendToNeighbourLoopBody(toRank,4,0,_vertexData.getCellDescriptionsIndex(4),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,h,level);
    #endif
  }
}

void exahype::Vertex::receiveNeighbourDataLoopBody(
    const int                                    fromRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    destCellDescriptionIndex,
    const bool                                   mergeWithReceivedData,
    const bool                                   receiveNeighbourMetadata,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level) {
  const auto src  = delineariseIndex2(srcScalar);
  const auto dest = delineariseIndex2(destScalar);
  assertion(tarch::la::countEqualEntries(dest, src) == (DIMENSIONS-1));

  if ( hasToReceiveMetadata(fromRank,srcScalar,destScalar,adjacentRanks) ) {
    solvers::Solver::BoundaryFaceInfo face(dest,src);
    const auto barycentre = computeFaceBarycentre(x,h,face._direction,dest);
    logDebug("receiveNeighbourDataLoopBody(...)","from rank="<<fromRank<<",barycentre="<<barycentre.toString()<<",level="<<level<<",adjacentRanks="<<adjacentRanks);

    exahype::MetadataHeap::HeapEntries receivedMetadata;
    receivedMetadata.clear();
    if ( receiveNeighbourMetadata ) {
      exahype::receiveNeighbourCommunicationMetadata(receivedMetadata,fromRank,barycentre,level);
      assertionEquals(receivedMetadata.size(),exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());
    }

    const bool validIndex = destCellDescriptionIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(destCellDescriptionIndex);
      for ( unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0; ) {
        auto* solver = solvers::RegisteredSolvers[solverNumber];
        const int begin = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
        const int end   = begin+exahype::NeighbourCommunicationMetadataPerSolver;

        switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          if ( receiveNeighbourMetadata ) {
            MetadataHeap::HeapEntries metadataPortion(receivedMetadata.begin()+begin,receivedMetadata.begin()+end);
            solvers::ADERDGSolver::mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,dest,src,barycentre);
          }
          if ( mergeWithReceivedData ) {
            static_cast<solvers::ADERDGSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
          } else {
            static_cast<solvers::ADERDGSolver*>(solver)->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
          }
          break;
        case solvers::Solver::Type::LimitingADERDG:
          if ( receiveNeighbourMetadata ) {
            MetadataHeap::HeapEntries metadataPortion(receivedMetadata.begin()+begin,receivedMetadata.begin()+end);
            solvers::ADERDGSolver::mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,dest,src,barycentre);
          }
          if ( mergeWithReceivedData ) {
            static_cast<solvers::LimitingADERDGSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
          } else {
            static_cast<solvers::LimitingADERDGSolver*>(solver)->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
          }
          break;
        case solvers::Solver::Type::FiniteVolumes:
          if ( mergeWithReceivedData ) {
            static_cast<solvers::FiniteVolumesSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,barycentre,level);
          } else {
            static_cast<solvers::FiniteVolumesSolver*>(solver)->dropNeighbourData(fromRank,barycentre,level);
          }
          break;
        default:
          assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          std::abort();
          break;
        }
      }
    }
  }
}

void exahype::Vertex::receiveNeighbourData(
    const int                                    fromRank,
    const bool                                   mergeWithReceivedData,
    const bool                                   isFirstIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    const bool receiveMetadata = !solvers::Solver::DisableMetadataExchangeDuringTimeSteps;

    const tarch::la::Vector<DIMENSIONS,int> lowerLeft(0);
    #if DIMENSIONS==3
    receiveNeighbourDataLoopBody(fromRank,4,0,_vertexData.getCellDescriptionsIndex(0),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
    receiveNeighbourDataLoopBody(fromRank,0,4,_vertexData.getCellDescriptionsIndex(4),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
    #endif
    receiveNeighbourDataLoopBody(fromRank,2,0,_vertexData.getCellDescriptionsIndex(0),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
    receiveNeighbourDataLoopBody(fromRank,0,2,_vertexData.getCellDescriptionsIndex(2),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
    receiveNeighbourDataLoopBody(fromRank,1,0,_vertexData.getCellDescriptionsIndex(0),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
    receiveNeighbourDataLoopBody(fromRank,0,1,_vertexData.getCellDescriptionsIndex(1),mergeWithReceivedData,receiveMetadata,getAdjacentRanks(),x,h,level);
  }
}
#endif
