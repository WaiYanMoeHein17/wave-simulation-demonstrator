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

// Definition of AMR-specific routines of the ADERDGSolver

#include "exahype/solvers/ADERDGSolver.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////

void exahype::solvers::ADERDGSolver::mergeWithEmptyNeighbourDuringMeshRefinement(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo, // corresponds to dest
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    const tarch::la::Vector<DIMENSIONS, int>&    posNeighbour,
    const bool                                   isAtBoundary,
    const tarch::la::Vector<DIMENSIONS, double>& barycentreFromVertex) {
  if ( isAtBoundary ) {
    mergeWithNeighbourMetadata(solverNumber,cellInfo,
          BoundaryStatus,BoundaryStatus,BoundaryStatus,
          pos,posNeighbour,barycentreFromVertex);
  } else {
    mergeWithNeighbourMetadata(solverNumber,cellInfo,
        EmptyStatus,EmptyStatus,EmptyStatus,
        pos,posNeighbour,barycentreFromVertex);
  }
}

void exahype::solvers::ADERDGSolver::mergeWithRefinementStatus(
    CellDescription& cellDescription,
    const int        faceIndex,
    const int        otherRefinementStatus) {
  cellDescription.setFacewiseRefinementStatus( faceIndex, otherRefinementStatus );
}

void
exahype::solvers::ADERDGSolver::updateCommunicationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  cellDescription.setCommunicationStatus(determineCommunicationStatus(cellDescription));
  assertion1(
      !isLeaf(cellDescription) || cellDescription.getCommunicationStatus()==LeafCommunicationStatus,
      cellDescription.toString());
}

int
exahype::solvers::ADERDGSolver::determineCommunicationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  if ( isLeaf(cellDescription) ) {
    return LeafCommunicationStatus;
  } else {
    int max = 0;
    for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      max = std::max( max, cellDescription.getFacewiseCommunicationStatus(i)-1 );
    }
    return max;
  }
}

void exahype::solvers::ADERDGSolver::mergeWithCommunicationStatus(
    CellDescription& cellDescription,
    const int faceIndex,
    const int otherCommunicationStatus) {
  assertion3(cellDescription.getCommunicationStatus()<=LeafCommunicationStatus,
             cellDescription.getCommunicationStatus(),otherCommunicationStatus,
             cellDescription.getCommunicationStatus());
  cellDescription.setFacewiseCommunicationStatus( faceIndex, otherCommunicationStatus );
}

void
exahype::solvers::ADERDGSolver::updateAugmentationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  cellDescription.setAugmentationStatus(determineAugmentationStatus(cellDescription));
  assertion1(
      !isParent(cellDescription) || cellDescription.getAugmentationStatus()==MaximumAugmentationStatus,
      cellDescription.toString());
}

int
exahype::solvers::ADERDGSolver::determineAugmentationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  if ( isParent(cellDescription) ) {
    return MaximumAugmentationStatus;
  } else {
    int max = 0;
    for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      max = std::max( max, cellDescription.getFacewiseAugmentationStatus(i)-1 );
    }
    return max;
  }
}

void exahype::solvers::ADERDGSolver::mergeWithAugmentationStatus(
    CellDescription& cellDescription,
    const int faceIndex,
    const int otherAugmentationStatus) {
  assertion3(
      cellDescription.getAugmentationStatus()<=MaximumAugmentationStatus,
      cellDescription.getAugmentationStatus(),otherAugmentationStatus,
      cellDescription.getAugmentationStatus());
  cellDescription.setFacewiseAugmentationStatus( faceIndex, otherAugmentationStatus );
}

void exahype::solvers::ADERDGSolver::updateRefinementStatus(CellDescription& cellDescription) const {
  if ( cellDescription.getLevel()==getMaximumAdaptiveMeshLevel() ) {
    int max = cellDescription.getRefinementFlag()  ?
        _refineOrKeepOnFineGrid : Erase;
    max = cellDescription.getRefinementStatus() == getMinRefinementStatusForTroubledCell() ?
        getMinRefinementStatusForTroubledCell() : max;

    for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      max = std::max( max, cellDescription.getFacewiseRefinementStatus(i)-1 );
    }
    cellDescription.setRefinementStatus(max);
  }
}


void exahype::solvers::ADERDGSolver::vetoParentCoarseningRequestIfNecessary(
    const CellDescription& fineGridCellDescription,
    CellDescription& coarseGridCellDescription) {
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool coarseGridCellDescriptionMighRequestCoarsening  =
      coarseGridCellDescription.getType()==CellDescription::Type::Parent ||
      coarseGridCellDescription.getType()==CellDescription::Type::ParentRequestsCoarseningA ||
      coarseGridCellDescription.getType()==CellDescription::Type::ParentRequestsCoarseningB;
  if (
      coarseGridCellDescriptionMighRequestCoarsening &&
      isLeaf(fineGridCellDescription)
  ) {
    if (
        fineGridCellDescription.getType()==CellDescription::Type::LeafInitiatesRefining ||
        fineGridCellDescription.getType()==CellDescription::Type::LeafRefines ||
        fineGridCellDescription.getRefinementStatus() >= Keep ||
        fineGridCellDescription.getPreviousRefinementStatus() >= Keep
    ) {
      coarseGridCellDescription.setType(CellDescription::Type::ParentChecked); // do not request coarsening again
    } else if (
        fineGridCellDescription.getRefinementStatus() == Pending ||
        fineGridCellDescription.getType()!=CellDescription::Type::LeafChecked ) {
      coarseGridCellDescription.setType(CellDescription::Type::Parent); // allows to request coarsening again
    }
  } else if (
      coarseGridCellDescriptionMighRequestCoarsening &&
      isParent(fineGridCellDescription)
  ) {
    if ( fineGridCellDescription.getType()==CellDescription::Type::ParentChecked ) { // do not request coarsening again
      coarseGridCellDescription.setType(CellDescription::Type::ParentChecked);
    } else {
      coarseGridCellDescription.setType(CellDescription::Type::Parent); // allows to request coarsening again
    }
  }
  lock.free();
}

// Allocate / deallocate

void exahype::solvers::ADERDGSolver::addNewCellDescription(
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const CellDescription::Type&                 cellType,
    const int                                    level,
    const int                                    parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const tarch::la::Vector<DIMENSIONS, double>& cellOffset) {
  assertion2(parentIndex == -1 || parentIndex != cellInfo._cellDescriptionsIndex, parentIndex, cellInfo._cellDescriptionsIndex);
  assertion2(parentIndex != cellInfo._cellDescriptionsIndex, parentIndex, cellInfo._cellDescriptionsIndex);
  logDebug("addNewCellDescription(...)","Add cell description: index="<<cellInfo._cellDescriptionsIndex<<",type="<<CellDescription::toString(cellType)<<",level="<<level<<",parentIndex="<<parentIndex << ",solver=" <<solverNumber);

  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),solverNumber,exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Background job completion monitoring (must be initialised with true)
  newCellDescription.setHasCompletedLastStep(true);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setParentIndex(parentIndex);
  newCellDescription.setParentType(static_cast<CellDescription::Type>(-1));
  newCellDescription.setLevel(level);

  newCellDescription.setAugmentationStatus(0);
  newCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
  newCellDescription.setCommunicationStatus(0);
  newCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
  if ( isLeaf(newCellDescription) ) {
    newCellDescription.setCommunicationStatus(LeafCommunicationStatus);
    newCellDescription.setFacewiseCommunicationStatus(LeafCommunicationStatus); // implicit conversion
    // TODO(Dominic): Make sure prolongation and restriction considers this.
  }
  newCellDescription.setNeighbourMergePerformed((signed char) 0/*implicit conversion*/);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Default field data indices
  newCellDescription.setSolutionIndex(-1);
  newCellDescription.setSolution(nullptr);
  newCellDescription.setPreviousSolutionIndex(-1);
  newCellDescription.setPreviousSolution(nullptr);
  newCellDescription.setUpdateIndex(-1);
  newCellDescription.setUpdate(nullptr);
  newCellDescription.setExtrapolatedPredictorIndex(-1);
  newCellDescription.setExtrapolatedPredictor(nullptr);
  newCellDescription.setFluctuationIndex(-1);
  newCellDescription.setFluctuation(nullptr);

  // Halo/Limiter meta data (oscillations identificator)
  newCellDescription.setRefinementFlag(false);
  newCellDescription.setRefinementStatus(Pending);
  newCellDescription.setPreviousRefinementStatus(Pending);
  newCellDescription.setFacewiseRefinementStatus(Pending);  // implicit conversion
  newCellDescription.setSolutionMinIndex(-1);
  newCellDescription.setSolutionMin(0);
  newCellDescription.setSolutionMaxIndex(-1);
  newCellDescription.setSolutionMax(0);

  // Compression
  newCellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);
  newCellDescription.setSolutionAveragesIndex(-1);
  newCellDescription.setSolutionAverages(nullptr);
  newCellDescription.setPreviousSolutionAveragesIndex(-1);
  newCellDescription.setPreviousSolutionAverages(nullptr);
  newCellDescription.setUpdateAveragesIndex(-1);
  newCellDescription.setUpdateAverages(nullptr);
  newCellDescription.setExtrapolatedPredictorAveragesIndex(-1);
  newCellDescription.setExtrapolatedPredictorAverages(nullptr);
  newCellDescription.setFluctuationAveragesIndex(-1);
  newCellDescription.setFluctuationAverages(nullptr);

  newCellDescription.setSolutionCompressedIndex(-1);
  newCellDescription.setSolutionCompressed(nullptr);
  newCellDescription.setPreviousSolutionAveragesIndex(-1);
  newCellDescription.setPreviousSolutionAverages(nullptr);
  newCellDescription.setUpdateCompressedIndex(-1);
  newCellDescription.setUpdateCompressed(nullptr);
  newCellDescription.setExtrapolatedPredictorCompressedIndex(-1);
  newCellDescription.setExtrapolatedPredictorCompressed(nullptr);
  newCellDescription.setFluctuationCompressedIndex(-1);
  newCellDescription.setFluctuationCompressed(nullptr);

  newCellDescription.setBytesPerDoFInExtrapolatedPredictor(-1);
  newCellDescription.setBytesPerDoFInFluctuation(-1);
  newCellDescription.setBytesPerDoFInPreviousSolution(-1);
  newCellDescription.setBytesPerDoFInSolution(-1);
  newCellDescription.setBytesPerDoFInUpdate(-1);

  newCellDescription.setCreation(CellDescription::Creation::NotSpecified);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  cellInfo._ADERDGCellDescriptions.push_back(newCellDescription);
  lock.free();
}

void exahype::solvers::ADERDGSolver::ensureNecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  // allocate solution
  if (
      holdsSolution(cellDescription) &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    // Allocate volume DoF
    const int dataPerCell = getDataPerCell(); // Only the solution and previousSolution store material parameters
    cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    cellDescription.setSolutionIndex        ( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionIndex(),"getPreviousSolutionIndex()");
    checkDataHeapIndex(cellDescription,cellDescription.getSolutionIndex(),"getSolutionIndex()");
    cellDescription.setPreviousSolution( getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).data() ) ;
    cellDescription.setSolution        ( getDataHeapEntries(cellDescription.getSolutionIndex()).data() );
    std::fill_n(static_cast<double*>(cellDescription.getPreviousSolution()),getDataPerCell(),std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getSolution()),getDataPerCell(),std::numeric_limits<double>::quiet_NaN());

    cellDescription.setSolutionCompressedIndex(-1);
    cellDescription.setSolutionCompressed(nullptr);
    cellDescription.setPreviousSolutionCompressedIndex(-1);
    cellDescription.setPreviousSolutionCompressed(nullptr);

    const int dataPerNode = getNumberOfVariables()+getNumberOfParameters();
    cellDescription.setPreviousSolutionAveragesIndex( DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );
    cellDescription.setSolutionAveragesIndex(         DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );
    checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionAveragesIndex(),"getPreviousSolutionAveragesIndex()");
    checkDataHeapIndex(cellDescription,cellDescription.getSolutionAveragesIndex(),"getSolutionAveragesIndex()");
    cellDescription.setPreviousSolutionAverages( getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex()).data() ) ;
    cellDescription.setSolutionAverages        ( getDataHeapEntries(cellDescription.getSolutionAveragesIndex()).data() ) ;
    std::fill_n(static_cast<double*>(cellDescription.getPreviousSolutionAverages()),dataPerNode,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getSolutionAverages()),dataPerNode,std::numeric_limits<double>::quiet_NaN());

    cellDescription.setCompressionState(CellDescription::Uncompressed);

    lock.free();
  }

  // allocate update and boundary arrays
  if (
      holdsFaceData(cellDescription) &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    // allocate update dof
    cellDescription.setUpdateIndex        ( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
    cellDescription.setUpdateAveragesIndex( DataHeap::getInstance().createData( getNumberOfVariables(), getNumberOfVariables() ) );
    cellDescription.setUpdateCompressedIndex(-1);
    cellDescription.setUpdateCompressed(nullptr);
    checkDataHeapIndex(cellDescription,cellDescription.getUpdateIndex(),"getUpdate()");
    checkDataHeapIndex(cellDescription,cellDescription.getUpdateAveragesIndex(),"getUpdateAverages()");
    cellDescription.setUpdate        ( getDataHeapEntries(cellDescription.getUpdateIndex()).data() ) ;
    cellDescription.setUpdateAverages( getDataHeapEntries(cellDescription.getUpdateAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getUpdate()),getUpdateSize(),std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getUpdateAverages()),getNumberOfVariables(),std::numeric_limits<double>::quiet_NaN());

    // extrapolated predictor
    const int dataPerBnd = getBndTotalSize();
    cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData(dataPerBnd, dataPerBnd) );
    cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
    cellDescription.setExtrapolatedPredictorCompressed(nullptr);
    const int boundaryData = (getNumberOfParameters()+getNumberOfVariables()) * DIMENSIONS_TIMES_TWO; //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setExtrapolatedPredictorAveragesIndex( DataHeap::getInstance().createData( boundaryData,  boundaryData  ) );
    checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedPredictorIndex(),"getExtrapolatedPredictor()");
    checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedPredictorAveragesIndex(),"getExtrapolatedPredictorAverages()");
    cellDescription.setExtrapolatedPredictor        ( getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).data() ) ;
    cellDescription.setExtrapolatedPredictorAverages( getDataHeapEntries(cellDescription.getExtrapolatedPredictorAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictor()),dataPerBnd,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages()),boundaryData,std::numeric_limits<double>::quiet_NaN());

    if (isUseViscousFlux()) {
      // gradients of extrapolated predictor
      const int gradientSizePerBnd = getBndGradQTotalSize();
      cellDescription.setExtrapolatedPredictorGradientIndex( DataHeap::getInstance().createData(gradientSizePerBnd, gradientSizePerBnd) );
      cellDescription.setExtrapolatedPredictorGradient( getDataHeapEntries(cellDescription.getExtrapolatedPredictorGradientIndex()).data() ) ;
      // touch the memory
      std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient()),gradientSizePerBnd,std::numeric_limits<double>::quiet_NaN());
    } else {
      cellDescription.setExtrapolatedPredictorGradientIndex(-1);
    }

    //
    // fluctuations
    //
    const int dofPerBnd  = getBndFluxTotalSize();
    cellDescription.setFluctuationIndex( DataHeap::getInstance().createData(dofPerBnd,  dofPerBnd) );
    cellDescription.setFluctuationCompressedIndex(-1);
    cellDescription.setFluctuationCompressed(nullptr);
    const int boundaryUnknowns = getNumberOfVariables() * DIMENSIONS_TIMES_TWO;     //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setFluctuationAveragesIndex( DataHeap::getInstance().createData( boundaryUnknowns, boundaryUnknowns ) );
    checkDataHeapIndex(cellDescription,cellDescription.getFluctuationIndex(),"getFluctuation()");
    checkDataHeapIndex(cellDescription,cellDescription.getFluctuationAveragesIndex(),"getFluctuationAverages()");
    cellDescription.setFluctuation        ( getDataHeapEntries(cellDescription.getFluctuationIndex()).data() ) ;
    cellDescription.setFluctuationAverages( getDataHeapEntries(cellDescription.getFluctuationAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getFluctuation()),dofPerBnd,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getFluctuationAverages()),boundaryUnknowns,std::numeric_limits<double>::quiet_NaN());

    // Allocate volume DoF for limiter (we need for every of the 2*DIMENSIONS faces an array of min values
    // and array of max values of the neighbour at this face).
    const int numberOfObservables = getDMPObservables();
    if ( numberOfObservables>0 ) {
      cellDescription.setSolutionMinIndex(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));
      cellDescription.setSolutionMaxIndex(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));
      checkDataHeapIndex(cellDescription,cellDescription.getSolutionMinIndex(),"getSolutionMinIndex()");
      checkDataHeapIndex(cellDescription,cellDescription.getSolutionMaxIndex(),"getSolutionMaxIndex()");
      cellDescription.setSolutionMin( getDataHeapEntries(cellDescription.getSolutionMinIndex()).data() ) ;
      cellDescription.setSolutionMax( getDataHeapEntries(cellDescription.getSolutionMaxIndex()).data() ) ;
      // touch the memory
      std::fill_n(static_cast<double*>(cellDescription.getSolutionMin()),numberOfObservables * DIMENSIONS_TIMES_TWO,std::numeric_limits<double>::quiet_NaN());
      std::fill_n(static_cast<double*>(cellDescription.getSolutionMax()),numberOfObservables * DIMENSIONS_TIMES_TWO,std::numeric_limits<double>::quiet_NaN());
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {

  if (
      !holdsSolution(cellDescription) &&
      DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())
  ) {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    if ( cellDescription.getSolutionIndex()>=0 ) {
      DataHeap::getInstance().deleteData(cellDescription.getSolutionIndex());
      assertion(cellDescription.getSolutionCompressedIndex()==-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getSolutionIndex()==-1);
      CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressedIndex());
    }

    DataHeap::getInstance().deleteData(cellDescription.getSolutionAveragesIndex());
    DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAveragesIndex());

    cellDescription.setPreviousSolutionIndex(-1);
    cellDescription.setPreviousSolution(nullptr);
    cellDescription.setSolutionIndex(-1);
    cellDescription.setSolution(nullptr);

    cellDescription.setPreviousSolutionAveragesIndex(-1);
    cellDescription.setPreviousSolutionAverages(nullptr);
    cellDescription.setSolutionAveragesIndex(-1);
    cellDescription.setSolutionAverages(nullptr);

    cellDescription.setPreviousSolutionCompressedIndex(-1);
    cellDescription.setPreviousSolutionCompressed(nullptr);
    cellDescription.setSolutionCompressedIndex(-1);
    cellDescription.setSolutionCompressed(nullptr);

    lock.free();
  }

  // deallocate update and boundary arrays
  if (
      !holdsFaceData(cellDescription) &&
      DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex())
  ) {
    // update
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()));
    if ( cellDescription.getUpdateIndex()>=0 ) {
      assertion(cellDescription.getUpdateCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getUpdateIndex());
      cellDescription.setUpdateIndex(-1);
      cellDescription.setUpdate(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getUpdateIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getUpdateCompressedIndex());
      cellDescription.setUpdateCompressedIndex(-1);
      cellDescription.setUpdateCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getUpdateAveragesIndex());
    cellDescription.setUpdateAveragesIndex(-1);
    cellDescription.setUpdateAverages(nullptr);

    // extrapolated predictor
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    if ( cellDescription.getExtrapolatedPredictorIndex()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
      assertion(cellDescription.getExtrapolatedPredictorCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorIndex());
      cellDescription.setExtrapolatedPredictorIndex(-1);
      cellDescription.setExtrapolatedPredictor(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getExtrapolatedPredictorIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorCompressedIndex());
      cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
      cellDescription.setExtrapolatedPredictorCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorAveragesIndex());
    cellDescription.setExtrapolatedPredictorAveragesIndex(-1);
    cellDescription.setExtrapolatedPredictorAverages(nullptr);

    // gradient of extrapolated predictor
    if (isUseViscousFlux() && cellDescription.getExtrapolatedPredictorGradientIndex() >= 0) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorGradientIndex()));

      DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorGradientIndex());
      cellDescription.setExtrapolatedPredictorGradientIndex(-1);
      cellDescription.setExtrapolatedPredictorGradient(nullptr);
    }

    // fluctuations
    if ( cellDescription.getFluctuationIndex()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));
      assertion(cellDescription.getFluctuationCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getFluctuationIndex());
      cellDescription.setFluctuationIndex(-1);
      cellDescription.setFluctuation(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getFluctuationIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getFluctuationCompressedIndex());
      cellDescription.setFluctuationCompressedIndex(-1);
      cellDescription.setFluctuationCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getFluctuationAveragesIndex());
    cellDescription.setFluctuationAveragesIndex(-1);
    cellDescription.setFluctuationAverages(nullptr);

    if ( getDMPObservables()>0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMinIndex()));
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMaxIndex()));
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMinIndex());
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMaxIndex());

      cellDescription.setSolutionMinIndex(-1);
      cellDescription.setSolutionMin(nullptr);
      cellDescription.setSolutionMaxIndex(-1);
      cellDescription.setSolutionMax(nullptr);
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::checkDataHeapIndex(const CellDescription& cellDescription, const int arrayIndex,const std::string arrayName) {
  assertion1(DataHeap::getInstance().isValidIndex(arrayIndex),cellDescription.toString());
  if ( arrayIndex < 0 ) {
    logError("checkDataHeapIndex(...)","The data heap array 'cellDescription."<<arrayName<<"' could not be allocated! Likely reason: Not enough memory available." <<
             " CellDescription="<<cellDescription.toString());
    std::abort();
  }
}

void exahype::solvers::ADERDGSolver::eraseCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion(Heap::getInstance().isValidIndex(cellDescriptionsIndex));
  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    auto *solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

    ADERDGSolver* aderdgSolver = nullptr;
    if (solver->getType()==Solver::Type::ADERDG) {
      aderdgSolver = static_cast<ADERDGSolver*>(solver);
    }
    else if (solver->getType()==Solver::Type::LimitingADERDG) {
      aderdgSolver =
          static_cast<LimitingADERDGSolver*>(solver)->getSolver().get();
    }
    assertion(aderdgSolver!=nullptr);

    p.setType(CellDescription::Type::Erased);
    aderdgSolver->ensureNoUnnecessaryMemoryIsAllocated(p);
  }

  Heap::getInstance().getData(cellDescriptionsIndex).clear();
}

// Operations


void exahype::solvers::ADERDGSolver::adjustSolutionDuringMeshRefinement(
    const int solverNumber,
    CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    synchroniseTimeStepping(cellDescription);

    const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
    if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
      cellDescription.setHasCompletedLastStep(false);
      peano::datatraversal::TaskSet( new AdjustSolutionDuringMeshRefinementJob(*this,cellDescription,isInitialMeshRefinement));
    } else {
      adjustSolutionDuringMeshRefinementBody(cellDescription,isInitialMeshRefinement);
    }
  }
}

void exahype::solvers::ADERDGSolver::adjustSolutionDuringMeshRefinementBody(
    CellDescription& cellDescription,
    const bool isInitialMeshRefinement) {
  #ifdef USE_ITAC
  VT_begin(adjustSolutionHandle);
  #endif

  if (
      cellDescription.getType()==CellDescription::Type::Leaf ||
      cellDescription.getType()==CellDescription::Type::LeafProlongates
  ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    if ( cellDescription.getType()==CellDescription::Type::LeafProlongates ) {
      prolongateVolumeData(cellDescription,isInitialMeshRefinement);
    }

    adjustSolution(cellDescription);
    markForRefinement(cellDescription);
    cellDescription.setType(CellDescription::Type::LeafChecked);
  }
  cellDescription.setHasCompletedLastStep(true);

  #ifdef USE_ITAC
  VT_end(adjustSolutionHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::updateStatusFlags(CellDescription& cellDescription) {
  const int preUpdateRefinementStatus    = cellDescription.getRefinementStatus();
  const int preUpdateAugmentationStatus  = cellDescription.getAugmentationStatus();
  const int preUpdateCommunicationStatus = cellDescription.getCommunicationStatus();

  updateCommunicationStatus(cellDescription);
  updateAugmentationStatus(cellDescription);
  updateRefinementStatus(cellDescription);

  if ( preUpdateRefinementStatus    != cellDescription.getRefinementStatus()    ) { AllSolversAreStable = false; }
  if ( preUpdateAugmentationStatus  != cellDescription.getAugmentationStatus()  ) { AllSolversAreStable = false; }
  if ( preUpdateCommunicationStatus != cellDescription.getCommunicationStatus() ) { AllSolversAreStable = false; }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int  solverNumber) {
  bool newComputeCell = false;

  // Fine grid cell based uniform mesh refinement.
  const int fineGridElement   = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
      fineGridElement==exahype::solvers::Solver::NotFound &&
      fineGridVerticesEnumerator.getLevel()==_coarsestMeshLevel
  ) {
    logDebug("progressMeshRefinementInEnterCell(...)","Add new uniform grid cell at centre="<<fineGridVerticesEnumerator.getCellCenter() <<", level="<<fineGridVerticesEnumerator.getLevel() << " for solver=" << solverNumber);

    addNewLeafCell(
        fineGridCell,fineGridVerticesEnumerator,
        mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex,
        solverNumber);
    newComputeCell = true;
  }
  else if ( fineGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
    ensureFineGridCoarseGridConsistency(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex()); // must come after refinement status update

    waitUntilCompletedLastStep(fineGridCellDescription,false,false); // wait for background jobs to complete

    // Update the status flags, allocate / deallocate memory
    updateStatusFlags(fineGridCellDescription);
    ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    finishRefinementOperation(fineGridCellDescription);
    decideOnRefinement(fineGridCellDescription);

    const int coarseGridElement = tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    if ( coarseGridElement != exahype::solvers::Solver::NotFound ) {
      CellDescription& coarseGridCellDescription = getCellDescription(
          coarseGridCell.getCellDescriptionsIndex(),coarseGridElement);
      vetoParentCoarseningRequestIfNecessary(fineGridCellDescription,coarseGridCellDescription);
    }

    // only assertions below
    // ...
    #ifdef Asserts
    const tarch::la::Vector<DIMENSIONS,double> center = fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize();
    #endif
    assertion5(Vertex::equalUpToRelativeTolerance(fineGridVerticesEnumerator.getCellCenter(),center),
        fineGridVerticesEnumerator.getCellCenter(),center,fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
    assertionEquals3(fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),fineGridVerticesEnumerator.getCellCenter(),fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize(),tarch::parallel::Node::getInstance().getRank());
  }

  // Coarse grid cell based adaptive mesh refinement operations.
  // Add new cells to the grid and veto erasing or erasing virtual children
  // requests if there are cells on the fine level.
  const int coarseGridElement = tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridElement);

    addNewVirtualCellIfVirtualRefiningRequested(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
        coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    newComputeCell |=
        addNewLeafIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
  }

  checkIfCellIsStable(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);

  return newComputeCell;
}

int exahype::solvers::ADERDGSolver::evaluateRefinementCriterion(
    const CellDescription& cellDescription, const double* const solution, const double& timeStamp) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf ||
             cellDescription.getType()==CellDescription::Type::LeafProlongates ||
             cellDescription.getType()==CellDescription::Type::ParentCoarsens
            ,cellDescription.toString());

  Solver::RefinementControl refinementControl =
      refinementCriterion(
          solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          timeStamp, // is done after the update
          cellDescription.getLevel());

  switch ( refinementControl ) {
    case exahype::solvers::Solver::RefinementControl::Refine:
        return _refineOrKeepOnFineGrid;
    case exahype::solvers::Solver::RefinementControl::Keep:
      return ( cellDescription.getLevel()==getMaximumAdaptiveMeshLevel() ) ? _refineOrKeepOnFineGrid : Keep;
    case exahype::solvers::Solver::RefinementControl::Erase:
      return Erase;
    default:
      logError("adjustSolutionDuringMeshRefinementBody(...)",
          "unknown refinement control value=" << static_cast<int>(refinementControl) <<
          ". Please check the return values of your refinement criterion.");
      std::abort();
      return Pending-1;
  }
}

void exahype::solvers::ADERDGSolver::markForRefinement(CellDescription& cellDescription) {
  const double* const solution = static_cast<double*>(cellDescription.getSolution());

  const int refinementStatus =
      evaluateRefinementCriterion(
          cellDescription,solution,
          cellDescription.getTimeStamp());

  if (
      refinementStatus==_refineOrKeepOnFineGrid &&
      cellDescription.getLevel()==getMaximumAdaptiveMeshLevel()
  ) {
    cellDescription.setRefinementFlag(true);
  }
  cellDescription.setRefinementStatus( std::max(cellDescription.getRefinementStatus(),refinementStatus) );
}

void exahype::solvers::ADERDGSolver::decideOnRefinement(CellDescription& fineGridCellDescription) {
  // top-down refining
  if (
     fineGridCellDescription.getType()==CellDescription::Type::LeafChecked &&
     fineGridCellDescription.getLevel()<getMaximumAdaptiveMeshLevel() &&
     fineGridCellDescription.getRefinementStatus() > Keep
  ) {
    fineGridCellDescription.setType(CellDescription::Type::LeafInitiatesRefining);
    AllSolversAreStable = false;
  }

  // bottom-up refining (halo refinement)
  if (
     fineGridCellDescription.getType()==CellDescription::Type::Virtual &&
     fineGridCellDescription.getLevel()==getMaximumAdaptiveMeshLevel() &&
     (fineGridCellDescription.getRefinementStatus() > Keep ||
     fineGridCellDescription.getPreviousRefinementStatus() > Keep)
  ) {
    Solver::SubcellPosition subcellPosition = amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(fineGridCellDescription);
    // top-most parent might have been erased
    if ( subcellPosition.parentElement != NotFound ) {
      CellDescription& topMostParent =
          getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
      tarch::multicore::Lock lock(ADERDGSolver::CoarseGridSemaphore);
      if ( topMostParent.getType()==CellDescription::Type::LeafChecked ) {
        topMostParent.setRefinementStatus(_refineOrKeepOnFineGrid);
        AllSolversAreStable = false;
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::addNewLeafCell(
    exahype::Cell& fineGridCell,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  logDebug("addNewCell(...)","Add new grid cell with center "<<fineGridVerticesEnumerator.getCellCenter() <<
              " at level "<<fineGridVerticesEnumerator.getLevel());

  AllSolversAreStable = false;

  CellInfo cellInfo = fineGridCell.addNewCellDescription(
      solverNumber,
      CellDescription::Type::Leaf,
      fineGridVerticesEnumerator.getLevel(),
      coarseGridCellDescriptionsIndex,
      fineGridVerticesEnumerator.getCellSize(),
      fineGridVerticesEnumerator.getVertexPosition());

  const int fineGridElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  CellDescription& fineGridCellDescription = cellInfo._ADERDGCellDescriptions[fineGridElement]; //TODO(Dominic): Multi-solvers: Might need to lock this?

  fineGridCellDescription.setRefinementStatus(Pending);
  fineGridCellDescription.setPreviousRefinementStatus(Keep);

  fineGridCellDescription.setCreation(CellDescription::Creation::UniformRefinement);
}

void exahype::solvers::ADERDGSolver::addNewVirtualCellIfVirtualRefiningRequested(
     exahype::Cell& fineGridCell,
     exahype::Vertex* const fineGridVertices,
     const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
     CellDescription& coarseGridCellDescription,
     const int coarseGridCellDescriptionsIndex) {
  // read from coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool addVirtualCell =
      (coarseGridCellDescription.getType()==CellDescription::Type::LeafChecked
      || isVirtual(coarseGridCellDescription)) &&
      coarseGridCellDescription.getAugmentationStatus() >= MinimumAugmentationStatusForVirtualRefining;
  lock.free();

  // work on fine grid
  const int fineGridElement = tryGetElement(
      fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());
  if (
      addVirtualCell &&
      fineGridElement == exahype::solvers::Solver::NotFound
  ) {
    AllSolversAreStable = false;

    fineGridCell.addNewCellDescription( // (EmptyDescendant),None
        coarseGridCellDescription.getSolverNumber(),
        CellDescription::Type::Virtual,
        fineGridVerticesEnumerator.getLevel(),
        coarseGridCellDescriptionsIndex,
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getVertexPosition());

    #ifdef Asserts
    const int fineGridElement = tryGetElement(fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());
    #endif
    assertion(fineGridElement!=exahype::solvers::Solver::NotFound);
  }
}

bool exahype::solvers::ADERDGSolver::addNewLeafIfRefinementRequested(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    CellDescription& coarseGridCellDescription,
    const int coarseGridCellDescriptionsIndex) {
  // read and modify coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  bool refiningOrRefiningRequested =
      coarseGridCellDescription.getType()==CellDescription::Type::LeafInitiatesRefining ||
      coarseGridCellDescription.getType()==CellDescription::Type::LeafRefines;
  if ( refiningOrRefiningRequested ) { // notifies that the fine grid cells exist
    coarseGridCellDescription.setType(CellDescription::Type::LeafRefines);
  }
  lock.free();

  // work on fine grid
  if ( refiningOrRefiningRequested ) {
    const int fineGridElement = tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());

    if ( fineGridElement==exahype::solvers::Solver::NotFound ) {
      addNewLeafCell(fineGridCell,fineGridVerticesEnumerator,
                 coarseGridCellDescriptionsIndex,
                 coarseGridCellDescription.getSolverNumber());
      CellDescription& fineGridCellDescription =
          getCellDescriptions(fineGridCell.getCellDescriptionsIndex()).back();
      fineGridCellDescription.setType(CellDescription::LeafProlongates);
      fineGridCellDescription.setCreation(CellDescription::Creation::AdaptiveRefinement);
    } else {
      CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);

      // only assertions
      #ifdef Parallel
      assertion4(fineGridCellDescription.getType()==CellDescription::Type::Virtual ||
                 fineGridCellDescription.getType()==CellDescription::Type::Leaf,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString(),
                 coarseGridCellDescriptionsIndex,
                 tarch::parallel::Node::getInstance().getRank());
      #else
      assertion2(fineGridCellDescription.getType()==CellDescription::Type::Virtual,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString());
      #endif
      assertion2(fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex,
                 fineGridCellDescription.toString(),coarseGridCellDescriptionsIndex);

      // the actual code begins here
      fineGridCellDescription.setType(CellDescription::Type::LeafProlongates);
      fineGridCellDescription.setCommunicationStatus(LeafCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      fineGridCellDescription.setRefinementStatus(Pending); // prevents that this cell is erased again
      fineGridCellDescription.setPreviousRefinementStatus(Keep);
      fineGridCellDescription.setCreation(CellDescription::Creation::AdaptiveRefinement);
    }
    return true;
  }
  return false;
}

void exahype::solvers::ADERDGSolver::prolongateVolumeData(
    CellDescription&       fineGridCellDescription,
    const bool initialGrid) {
  const int coarseGridElement =
      tryGetElement(fineGridCellDescription.getParentIndex(),fineGridCellDescription.getSolverNumber());
  assertion1(coarseGridElement!=exahype::solvers::Solver::NotFound,fineGridCellDescription.toString());
  CellDescription& coarseGridCellDescription =
      getCellDescription(fineGridCellDescription.getParentIndex(),coarseGridElement);

  tarch::la::Vector<DIMENSIONS,int> subcellIndex =
      amr::computeSubcellIndex(
          fineGridCellDescription.getOffset(),
          fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

  const int levelFine = fineGridCellDescription.getLevel();
  const int levelCoarse = coarseGridCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  // current solution
  double* solutionFine   = static_cast<double*>(fineGridCellDescription.getSolution());
  double* solutionCoarse = static_cast<double*>(coarseGridCellDescription.getSolution());
  volumeUnknownsProlongation(
      solutionFine,solutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  // previous solution
  assertion(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolutionIndex()));
  double* previousSolutionFine   = static_cast<double*>(fineGridCellDescription.getPreviousSolution());
  double* previousSolutionCoarse = static_cast<double*>(coarseGridCellDescription.getPreviousSolution());
  volumeUnknownsProlongation(
      previousSolutionFine,previousSolutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  fineGridCellDescription.setTimeStamp(coarseGridCellDescription.getTimeStamp());
  fineGridCellDescription.setTimeStepSize(coarseGridCellDescription.getTimeStepSize());
}

// TODO(Dominic): Add to docu: Do not compute gradient via cell-centered value as virtual cell might be surrounded by empty cells
bool exahype::solvers::ADERDGSolver::checkIfStatusFlaggingHasConverged(
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& neighbourFlags,
    const int                                          minValueToConsider) {
  bool flaggingHasConverged = true;
  for (int d=0; d<DIMENSIONS; d++) {
    if ( neighbourFlags(2*d+1) >= minValueToConsider && neighbourFlags(2*d+0) >= minValueToConsider ) {
      flaggingHasConverged &=
          std::abs(neighbourFlags(2*d+1) - neighbourFlags(2*d+0)) <= 2;
    }
  }
  return flaggingHasConverged;
}

void exahype::solvers::ADERDGSolver::checkIfCellIsStable(
        exahype::Cell&                       fineGridCell,
        exahype::Vertex* const               fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        const int                            solverNumber) const {
  const int element = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( element!=exahype::solvers::Solver::NotFound ) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);
    // compute flagging gradients in inside cells
    bool flaggingHasConverged =
        checkIfStatusFlaggingHasConverged(cellDescription.getFacewiseAugmentationStatus(),0)  &&
        checkIfStatusFlaggingHasConverged(cellDescription.getFacewiseCommunicationStatus(),0);
    // refinement status is only spread on finest level
    if (
        ( isLeaf(cellDescription) || isVirtual(cellDescription) ) &&
        cellDescription.getLevel() == getMaximumAdaptiveMeshLevel()
    ) {
      flaggingHasConverged &= checkIfStatusFlaggingHasConverged(cellDescription.getFacewiseRefinementStatus(),Erase);
    }

    const bool noCoarseGridLeafCellMustRequireRefinement =
        !(isLeaf(cellDescription) &&
        cellDescription.getLevel() < getMaximumAdaptiveMeshLevel() &&
        cellDescription.getRefinementStatus()>=_refineOrKeepOnFineGrid);
    const bool noFinestGridVirtualCellMustBeReplacedByLeafCell =
        !(isVirtual(cellDescription) &&
        cellDescription.getLevel() == getMaximumAdaptiveMeshLevel() &&
        cellDescription.getRefinementStatus() > Keep);
    const bool noLeafCellIsNextToAnEmptyCell =
        !(isLeaf(cellDescription) &&
        tarch::la::min(cellDescription.getFacewiseRefinementStatus())==ADERDGSolver::EmptyStatus);

    const bool stable =
      flaggingHasConverged
      &&
      (cellDescription.getType()==CellDescription::Type::LeafChecked  ||
      cellDescription.getType()==CellDescription::Type::ParentChecked ||
      cellDescription.getType()==CellDescription::Type::Virtual)
      &&
      noLeafCellIsNextToAnEmptyCell
      &&
      noCoarseGridLeafCellMustRequireRefinement
      &&
      noFinestGridVirtualCellMustBeReplacedByLeafCell;

    if ( !stable ) {
      AllSolversAreStable = false;

      #ifdef MonitorMeshRefinement
      logInfo("checkIfCellIsStable(...)","cell has not attained stable state (yet):");
      logInfo("checkIfCellIsStable(...)","type="<<cellDescription.toString(cellDescription.getType()));
      logInfo("checkIfCellIsStable(...)","offset="<<cellDescription.getOffset());
      logInfo("checkIfCellIsStable(...)","size="<<cellDescription.getSize());
      logInfo("checkIfCellIsStable(...)","level="<<cellDescription.getLevel());
      logInfo("checkIfCellIsStable(...)","flaggingHasConverged="<<flaggingHasConverged);
      logInfo("checkIfCellIsStable(...)","noLeafCellIsNextToAnEmptyCell="<<noLeafCellIsNextToAnEmptyCell);
      logInfo("checkIfCellIsStable(...)","noCoarseGridLeafCellMustRequireRefinement="<<noCoarseGridLeafCellMustRequireRefinement);
      logInfo("checkIfCellIsStable(...)","noFinestGridVirtualCellMustBeReplacedByLeafCell="<<noFinestGridVirtualCellMustBeReplacedByLeafCell);
      logInfo("checkIfCellIsStable(...)","getAugmentationStatus="<<cellDescription.getAugmentationStatus());
      logInfo("checkIfCellIsStable(...)","getCommunicationStatus="<<cellDescription.getCommunicationStatus());
      logInfo("checkIfCellIsStable(...)","getRefinementStatus="<<cellDescription.getRefinementStatus());
      logInfo("checkIfCellIsStable(...)","getFacewiseAugmentationStatus="<<cellDescription.getFacewiseAugmentationStatus());
      logInfo("checkIfCellIsStable(...)","getFacewiseCommunicationStatus="<<cellDescription.getFacewiseCommunicationStatus());
      logInfo("checkIfCellIsStable(...)","getFacewiseRefinementStatus="<<cellDescription.getFacewiseRefinementStatus());
      logInfo("checkIfCellIsStable(...)","cellDescription.getCreation="<<cellDescription.toString(cellDescription.getCreation()));
      logInfo("checkIfCellIsStable(...)","fineGridCell.getCellDescriptionsIndex="<<fineGridCell.getCellDescriptionsIndex());
      logInfo("checkIfCellIsStable(...)","solver.getCoarsestMeshLevel="<<getCoarsestMeshLevel());
      logInfo("checkIfCellIsStable(...)","solver.getMaximumAdaptiveMeshLevel="<<getMaximumAdaptiveMeshLevel());
      logInfo("checkIfCellIsStable(...)","solver.getMaximumAdaptiveMeshDepth="<<getMaximumAdaptiveMeshDepth());
      #endif
    }
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  const int fineGridElement = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( fineGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);

    // start or finish collective operations
    progressCoarseningOperationsInLeaveCell(fineGridCellDescription);

    // skip remainder if the refinement criterion has not been evaluated yet for a Cell
    // Reading the refinement request might result into data race but this is accepted at this point
    // as we only read and not write
    CellDescription::Type parentType = CellDescription::Type::Erased;
    int parentAugmentationStatus     = 0;
    const int coarseGridElement = tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
      assertion3(fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),fineGridCellDescription.toString(),fineGridCell.toString(),coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.

      CellDescription& coarseGridCellDescription = getCellDescription(
          fineGridCellDescription.getParentIndex(),coarseGridElement);
      assertion1(fineGridCellDescription.getSolverNumber()==coarseGridCellDescription.getSolverNumber(),fineGridCellDescription.toString());

      restrictVolumeDataIfParentCoarsens(
          fineGridCellDescription,coarseGridCellDescription);

      parentType               = coarseGridCellDescription.getType();
      parentAugmentationStatus = coarseGridCellDescription.getAugmentationStatus();
    }

    eraseCellDescriptionIfNecessary(
        fineGridCellDescription,
        fineGridCell.getCellDescriptionsIndex(),
        fineGridElement,
        parentType,
        parentAugmentationStatus);
  }

  checkIfCellIsStable(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::ADERDGSolver::eraseOrRefineAdjacentVertices(
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
    const int isValidIndex =
        cellDescriptionsIndex > 0 &&
        (!checkThoroughly ||
        Heap::getInstance().getInstance().isValidIndex(cellDescriptionsIndex));

    if ( isValidIndex ) {
      const int element = tryGetElement(cellDescriptionsIndex,solverNumber);
      if ( element!=NotFound ) {
        CellDescription& cellDescription = getCellDescription(
            cellDescriptionsIndex,element);

        checkSuccessful = Vertex::equalUpToRelativeTolerance(cellDescription.getOffset(), cellOffset);
        if (
            !checkThoroughly || checkSuccessful
        ) {
          bool refineAdjacentVertex =
              isParent(cellDescription) ||
              cellDescription.getType()==CellDescription::Type::LeafInitiatesRefining ||
              cellDescription.getType()==CellDescription::Type::LeafRefines ||
              cellDescription.getAugmentationStatus() >= MinimumAugmentationStatusForRefining ||
              cellDescription.getRefinementStatus() > Keep ||
              cellDescription.getPreviousRefinementStatus() > Keep;
          refineAdjacentVertex &= level < getMaximumAdaptiveMeshLevel();

          bool eraseAdjacentVertex =
              (cellDescription.getType()==CellDescription::Type::LeafChecked ||
              cellDescription.getType()==CellDescription::Type::Virtual)
              &&
              cellDescription.getAugmentationStatus() < MinimumAugmentationStatusForRefining // TODO(Dominic): Probably can tune here. This is chosen to large
              &&
              cellDescription.getRefinementStatus() <= Keep &&
              cellDescription.getPreviousRefinementStatus() <= Keep;

          if (refineAdjacentVertex) {
            return RefinementControl::Refine;
          } else if (eraseAdjacentVertex) {
            return RefinementControl::Erase;
          } else {
            return RefinementControl::Keep;
          }
        } else {
          return RefinementControl::Keep; // ?
        }
      } else {
        return RefinementControl::Erase;
      }
    } else {
      return RefinementControl::Keep;
    }
  }
}

void exahype::solvers::ADERDGSolver::prepareVolumeDataRestriction(
    CellDescription& cellDescription) const {
  double* solution = static_cast<double*>(cellDescription.getSolution());
  double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
  std::fill_n(solution,getDataPerCell(),0.0);
  std::fill_n(previousSolution,getDataPerCell(),0.0);

  const int DMPObservables = getDMPObservables() > 0;
  if ( DMPObservables > 0 ) {
    double* DMPObservablesMin = static_cast<double*>(cellDescription.getSolutionMin());
    double* DMPObservablesMax = static_cast<double*>(cellDescription.getSolutionMax());
    std::fill_n(DMPObservablesMin,DIMENSIONS_TIMES_TWO*DMPObservables,+std::numeric_limits<double>::infinity());
    std::fill_n(DMPObservablesMax,DIMENSIONS_TIMES_TWO*DMPObservables,-std::numeric_limits<double>::infinity());
  }
}

void exahype::solvers::ADERDGSolver::changeLeafToParent(CellDescription& cellDescription) {
  assertion1(
      cellDescription.getType()==CellDescription::Type::LeafRefines ||
      cellDescription.getType()==CellDescription::Type::ParentCoarsens
      ,cellDescription.toString());
  cellDescription.setType(CellDescription::Type::ParentChecked);
  cellDescription.setAugmentationStatus(MaximumAugmentationStatus);
  cellDescription.setRefinementStatus(Keep);
  cellDescription.setPreviousRefinementStatus(Keep);
  cellDescription.setCommunicationStatus(0);
  cellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
  cellDescription.setFacewiseRefinementStatus(Pending);
  cellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
}

void exahype::solvers::ADERDGSolver::finishRefinementOperation(
     CellDescription& fineGridCellDescription) {
  if ( fineGridCellDescription.getType() == CellDescription::Type::LeafRefines ) {
    changeLeafToParent(fineGridCellDescription);
  }
}

bool exahype::solvers::ADERDGSolver::markPreviousParentForRefinement(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::ParentCoarsens,cellDescription.toString());
    double* solution = static_cast<double*>(cellDescription.getSolution());
    adjustSolution(solution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getTimeStamp(),
          cellDescription.getTimeStepSize());

    double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
    adjustSolution(previousSolution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getPreviousTimeStamp(),
          cellDescription.getPreviousTimeStepSize());

    cellDescription.setRefinementStatus(
      evaluateRefinementCriterion(
        cellDescription,solution,cellDescription.getTimeStamp())
    );
    cellDescription.setPreviousRefinementStatus(
        evaluateRefinementCriterion(
            cellDescription,previousSolution,cellDescription.getPreviousTimeStamp())
    );

    return cellDescription.getRefinementStatus()        <= Keep &&
           cellDescription.getPreviousRefinementStatus()<= Keep;
}

void exahype::solvers::ADERDGSolver::progressCoarseningOperationsInLeaveCell(CellDescription& fineGridCellDescription) {
  switch ( fineGridCellDescription.getType() ) {
    case CellDescription::Type::Parent:
      fineGridCellDescription.setType(CellDescription::Type::ParentRequestsCoarseningA);
      break; // in leaveCell since we might spawn refinement crit. eval as background job
    case CellDescription::Type::ParentRequestsCoarseningA:
      fineGridCellDescription.setType(CellDescription::Type::ParentRequestsCoarseningB);
      break; // in leaveCell since we might spawn refinement crit. eval as background job
    case CellDescription::Type::ParentRequestsCoarseningB:
      // requests was not vetoed
      fineGridCellDescription.setType(CellDescription::Type::ParentCoarsens);
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      fineGridCellDescription.setAugmentationStatus(0);
      fineGridCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
      fineGridCellDescription.setCommunicationStatus(LeafCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(LeafCommunicationStatus); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      prepareVolumeDataRestriction(fineGridCellDescription);
      break;
    case CellDescription::Type::ParentCoarsens:
      if ( markPreviousParentForRefinement(fineGridCellDescription) ) { // evaluate refinement criterion now that fine grid cells have restricted their data
        fineGridCellDescription.setType(CellDescription::Type::LeafChecked);
        fineGridCellDescription.setRefinementStatus(Keep);
      } else {
        changeLeafToParent(fineGridCellDescription);
        assertion1(fineGridCellDescription.getType()==CellDescription::Type::ParentChecked,fineGridCellDescription.toString());
      }
      break;
    default:
      break;
  }
}

// TODO(Dominic): Assumes existence of parent
// Must not exist in case of

void exahype::solvers::ADERDGSolver::eraseCellDescriptionIfNecessary(
    CellDescription&             cellDescription,
    const int                    cellDescriptionsIndex,
    const int                    fineGridElement,
    const CellDescription::Type& parentType,
    const int                    parentAugmentationStatus) const {
  bool eraseVirtualCell =
      cellDescription.getType() == CellDescription::Type::Virtual &&
      parentType==CellDescription::Type::Erased;
  eraseVirtualCell |=
      parentAugmentationStatus<MinimumAugmentationStatusForVirtualRefining
      &&
      cellDescription.getType() == CellDescription::Type::Virtual                                    &&
      cellDescription.getCommunicationStatus() < MinimumCommunicationStatusForNeighbourCommunication
      &&
      (cellDescription.getLevel() != getMaximumAdaptiveMeshLevel() ||
      cellDescription.getRefinementStatus() < Keep);

  const bool erasingLeafCell = // see progressCoarseningOperationsInLeaveCell
      cellDescription.getType() == CellDescription::Type::LeafChecked &&
      parentType == CellDescription::Type::LeafChecked;

  if ( erasingLeafCell || eraseVirtualCell ) {
    AllSolversAreStable = false;

    cellDescription.setType(CellDescription::Erased);
    cellDescription.setCommunicationStatus(0);
    cellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridElement);
  }
}

void exahype::solvers::ADERDGSolver::restrictVolumeDataIfParentCoarsens(
    const CellDescription& fineGridCellDescription,
    const CellDescription& coarseGridCellDescription) {
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool restrictVolumeData =
      coarseGridCellDescription.getType()==CellDescription::Type::ParentCoarsens;
  lock.free();

  if ( restrictVolumeData ) {
    tarch::la::Vector<DIMENSIONS,int> subcellIndex =
        exahype::amr::computeSubcellIndex(
            fineGridCellDescription.getOffset(),
            fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

    // restrict values.
    assertion1(fineGridCellDescription.getRefinementStatus()==Erase,fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getSolutionIndex()),fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(coarseGridCellDescription.getSolutionIndex()),coarseGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolutionIndex()),fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(coarseGridCellDescription.getPreviousSolutionIndex()),coarseGridCellDescription.toString());

    const int levelFine   = fineGridCellDescription.getLevel();
    const int levelCoarse = coarseGridCellDescription.getLevel();
    assertion(levelCoarse < levelFine);

    if ( !DataHeap::getInstance().isValidIndex(fineGridCellDescription.getSolutionIndex()) ) {
      logError("restrictVolumeData(..)","solution not valid for cell="<<fineGridCellDescription.toString());
      std::abort();
    }

    // restrict current solution
    tarch::multicore::Lock lock(RestrictionSemaphore);
    double* solutionFine   = static_cast<double*>(fineGridCellDescription.getSolution());
    double* solutionCoarse = static_cast<double*>(coarseGridCellDescription.getSolution());
    volumeUnknownsRestriction(
        solutionCoarse,solutionFine,
        levelCoarse,levelFine,
        subcellIndex);

    // restrict next solution
    double* previousSolutionFine   = static_cast<double*>(fineGridCellDescription.getPreviousSolution());
    double* previousSolutionCoarse = static_cast<double*>(coarseGridCellDescription.getPreviousSolution());
    volumeUnknownsRestriction(
        previousSolutionCoarse,previousSolutionFine,
        levelCoarse,levelFine,
        subcellIndex);

    // Restrict the observables min and max
    const int numberOfDMPObservables = getDMPObservables();
    if ( numberOfDMPObservables > 0 ) {
      double* coarseGridObservablesMin = static_cast<double*>(coarseGridCellDescription.getSolutionMin());
      double* coarseGridObservablesMax = static_cast<double*>(coarseGridCellDescription.getSolutionMax());
      const double* const fineGridObservablesMin = static_cast<double*>(fineGridCellDescription.getSolutionMin());
      const double* const fineGridObservablesMax = static_cast<double*>(fineGridCellDescription.getSolutionMax());

      for ( int i = 0; i < DIMENSIONS_TIMES_TWO*numberOfDMPObservables; i++ ) {
        coarseGridObservablesMin[i] = std::min( coarseGridObservablesMin[i], fineGridObservablesMin[i] );
        coarseGridObservablesMax[i] = std::max( coarseGridObservablesMax[i], fineGridObservablesMax[i] );
      }
    }

    // TODO(Dominic): What to do with the time step data for anarchic time stepping?
    // Tobias proposed some waiting procedure. Until they all have reached
    // a certain time level.
    //  coarseGridCellDescription.setTimeStamp(fineGridCellDescription.getTimeStamp());
    //  coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    //  coarseGridCellDescription.setTimeStepSize(fineGridCellDescription.getTimeStepSize());
    //  coarseGridCellDescription.setPredictorTimeStepSize(fineGridCellDescription.getPredictorTimeStepSize());
    // TODO(Dominic): Reconsider for anarchic time stepping.
    // coarseGridCellDescription.setTimeStamp(fineGridCellDescription.getTimeStamp());
    // coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::ensureFineGridCoarseGridConsistency(
    CellDescription& fineGridCellDescription,
    const int coarseGridCellDescriptionsIndex) {
  fineGridCellDescription.setParentIndex(mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex);
  if ( isValidCellDescriptionIndex(coarseGridCellDescriptionsIndex) ) {
    fineGridCellDescription.setParentIndex(coarseGridCellDescriptionsIndex);
  }
}

void exahype::solvers::ADERDGSolver::finaliseStateUpdates(
      const int solverNumber,
      CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=exahype::solvers::Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    assertion1(
        cellDescription.getType()==CellDescription::Type::LeafChecked   ||
        cellDescription.getType()==CellDescription::Type::ParentChecked ||
        cellDescription.getType()==CellDescription::Type::Virtual,
        cellDescription.toString());

    if ( cellDescription.getType()==CellDescription::Type::LeafChecked ) {
      cellDescription.setType(CellDescription::Type::Leaf);
    } else if ( cellDescription.getType()==CellDescription::Type::ParentChecked ) {
      cellDescription.setType(CellDescription::Type::Parent);
    }

    assertion1(
        cellDescription.getType()==CellDescription::Type::Leaf          ||
        cellDescription.getType()==CellDescription::Type::Parent        ||
        cellDescription.getType()==CellDescription::Type::Virtual,
        cellDescription.toString());

    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      validateCellDescriptionData(cellDescription,cellDescription.getTimeStamp()>0,false,true,"finaliseStateUpdates");
    }
    cellDescription.setRefinementFlag(false);

    if ( getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested ) {
      cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
    }
  }
}

//////////////////////////
/// Forks and Joins
//////////////////////////


#ifdef Parallel

const int exahype::solvers::ADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 2;

/**
 * After the forking the master's cell descriptions
 * are not accessed by enterCell(...) on the master
 * anymore. However we still use them as buffer for saving data.
 */
bool exahype::solvers::ADERDGSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const bool                                    fromWorkerSide,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    logDebug("sendCellDescriptions(...)","send "<< Heap::getInstance().getData(cellDescriptionsIndex).size()<<
        " cell descriptions to rank "<<toRank<<" (x="<< x.toString() << ",level="<< level << ")");
    bool oneSolverRequiresVerticalCommunication = false;
    for (auto& cellDescription : Heap::getInstance().getData(cellDescriptionsIndex)) {
      // wait for background jobs to complete
      if ( !cellDescription.getHasCompletedLastStep() ) {
        peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
      }
      int numberOfBackgroundJobsToProcess = 1;
      while ( !cellDescription.getHasCompletedLastStep() ) {
        tarch::multicore::jobs::processBackgroundJobs(numberOfBackgroundJobsToProcess);
        numberOfBackgroundJobsToProcess++;
      }
      oneSolverRequiresVerticalCommunication &=
          cellDescription.getType()==CellDescription::Type::Virtual &&
          cellDescription.getAugmentationStatus() >= MinimumAugmentationStatusForVirtualRefining;
    }
    Heap::getInstance().sendData(cellDescriptionsIndex,toRank,x,level,messageType);
    return oneSolverRequiresVerticalCommunication;
  }
  else {
    logDebug("sendCellDescriptions(...)","send "
        " empty cell descriptions to rank "<<toRank<<" (x="<< x.toString() << ",level="<< level << ")");

    sendEmptyCellDescriptions(toRank,messageType,x,level);
    return false;
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  logDebug("sendEmptyCellDescriptions(...)","send empty message to " <<
          "rank "<<toRank <<
          " at (center="<< x.toString() <<
          ",level="<< level << ")");

  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::receiveCellDescriptions(
    const int                                    fromRank,
    exahype::Cell&                               localCell,
    const peano::heap::MessageType&              messageType,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  Heap::getInstance().receiveData(
      localCell.getCellDescriptionsIndex(),fromRank,x,level,messageType);

  logDebug("mergeCellDescriptionsWithRemoteData(...)","received " <<
          Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).size() <<
          " cell descriptions for cell (centre="<< x.toString() << ", level="<< level << ")");

  for (auto& cellDescription : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
    resetIndicesAndFlagsOfReceivedCellDescription(
        cellDescription,mappings::LevelwiseAdjacencyBookkeeping::InvalidAdjacencyIndex);
  }
}

void exahype::solvers::ADERDGSolver::resetIndicesAndFlagsOfReceivedCellDescription(
    CellDescription& cellDescription,const int parentIndex) {
  cellDescription.setParentIndex(parentIndex);

  // Default field data indices
  cellDescription.setSolutionIndex(-1);
  cellDescription.setSolution(nullptr);
  cellDescription.setPreviousSolutionIndex(-1);
  cellDescription.setPreviousSolution(nullptr);
  cellDescription.setUpdateIndex(-1);
  cellDescription.setUpdate(nullptr);
  cellDescription.setExtrapolatedPredictorIndex(-1);
  cellDescription.setExtrapolatedPredictor(nullptr);
  cellDescription.setFluctuationIndex(-1);
  cellDescription.setFluctuation(nullptr);

  // Limiter meta data (oscillations identificator)
  cellDescription.setSolutionMinIndex(-1);
  cellDescription.setSolutionMin(nullptr);
  cellDescription.setSolutionMaxIndex(-1);
  cellDescription.setSolutionMax(nullptr);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
  cellDescription.setExtrapolatedPredictorCompressed(nullptr);
  cellDescription.setFluctuationCompressedIndex(-1);
  cellDescription.setFluctuationCompressed(nullptr);
  cellDescription.setSolutionCompressedIndex(-1);
  cellDescription.setSolutionCompressed(nullptr);
  cellDescription.setPreviousSolutionCompressedIndex(-1);
  cellDescription.setPreviousSolutionCompressed(nullptr);
  cellDescription.setUpdateCompressedIndex(-1);
  cellDescription.setUpdateCompressed(nullptr);

  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setUpdateAveragesIndex(-1);
  cellDescription.setUpdateAverages(nullptr);
  cellDescription.setExtrapolatedPredictorAveragesIndex(-1);
  cellDescription.setExtrapolatedPredictorAverages(nullptr);
  cellDescription.setFluctuationAveragesIndex(-1);
  cellDescription.setFluctuationAverages(nullptr);
  cellDescription.setBytesPerDoFInPreviousSolution(-1);
  cellDescription.setBytesPerDoFInSolution(-1);
  cellDescription.setBytesPerDoFInUpdate(-1);
  cellDescription.setBytesPerDoFInExtrapolatedPredictor(-1);
  cellDescription.setBytesPerDoFInFluctuation(-1);

  // reset the facewise flags
  cellDescription.setFacewiseAugmentationStatus(0);
  cellDescription.setFacewiseCommunicationStatus(0);

  cellDescription.setFacewiseRefinementStatus(Pending);

  cellDescription.setCreation(CellDescription::Creation::ReceivedDueToForkOrJoin);

  // background jobs
  cellDescription.setHasCompletedLastStep(true);
}

void exahype::solvers::ADERDGSolver::ensureOnlyNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  auto* solver = RegisteredSolvers[cellDescription.getSolverNumber()];
  switch (solver->getType()) {
  case exahype::solvers::Solver::Type::ADERDG:
    static_cast<ADERDGSolver*>(solver)->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<ADERDGSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::LimitingADERDG:
    static_cast<LimitingADERDGSolver*>(solver)->
    getSolver()->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<LimitingADERDGSolver*>(solver)->
        getSolver()->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::FiniteVolumes:
    assertionMsg(false,"Solver type not supported!");
    break;
  }
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::ADERDGSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

////////////////////////////////////
// MASTER <=> WORKER
////////////////////////////////////

void exahype::solvers::ADERDGSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  // coarse grid based operations
  const int coarseGridCellDescriptionsIndex = coarseGridCell.getCellDescriptionsIndex();
  const int coarseGridCellElement = tryGetElement(coarseGridCellDescriptionsIndex,solverNumber);

  if ( coarseGridCellElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement);

    addNewVirtualCellIfVirtualRefiningRequested(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
        coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    bool addedNewCell =
        addNewLeafIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

    if ( element!=exahype::solvers::Solver::NotFound ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      // !!! Use the parent type argument to communicate the coarse grid cell type to the worker rank
      fineGridCellDescription.setParentType(coarseGridCellDescription.getType());

      if ( isVirtual(fineGridCellDescription) ) {
        Solver::SubcellPosition subcellPosition =  amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(fineGridCellDescription);
        CellDescription& topMostParentCellDescription =
            getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
        if ( isLeaf(topMostParentCellDescription) ) {
           logDebug( "progressMeshRefinementInPrepareSendToWorker(...)"," try to refine parent " << topMostParentCellDescription.toString());
           topMostParentCellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
        }
      }
    }

    if ( addedNewCell ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      prolongateVolumeData(fineGridCellDescription,getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested);
      assertion1( fineGridCellDescription.getType()==CellDescription::Type::LeafProlongates,fineGridCellDescription.toString());
    }
  }
}

void exahype::solvers::ADERDGSolver::sendDataToWorkerIfProlongating(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);

  // send out the data
  if ( fineGridCellDescription.getType()==CellDescription::Type::LeafProlongates ) {
    logDebug( "sendDataToWorkerIfProlongating(...)","send prolongated solution to rank "<<workerRank<< " at x="<<x.toString()<< ",level="<<level << " cell="<<fineGridCellDescription.toString());

    sendDataToWorkerOrMasterDueToForkOrJoin(workerRank,cellDescriptionsIndex,element,
        peano::heap::MessageType::MasterWorkerCommunication,x,level);
  }
}

void exahype::solvers::ADERDGSolver::receiveDataFromMasterIfProlongating(
  const int masterRank,
    const int receivedCellDescriptionsIndex,
  const int receivedElement,
  const tarch::la::Vector<DIMENSIONS,double>& x,
  const int level) const {
  CellDescription& receivedCellDescription = getCellDescription(receivedCellDescriptionsIndex,receivedElement);

  if ( receivedCellDescription.getType()==CellDescription::Type::LeafProlongates ) {
    logDebug( "receiveDataFromMasterIfProlongating(...)","receiving prolongated solution from rank "<<masterRank<< " at x="<<x.toString()<< ",level="<<level);

    mergeWithWorkerOrMasterDataDueToForkOrJoin(
        masterRank,receivedCellDescriptionsIndex,receivedElement,
        peano::heap::MessageType::MasterWorkerCommunication,x,level);
  }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
  auto& receivedCellDescriptions = getCellDescriptions(receivedCellDescriptionsIndex);
  assertion1( isValidCellDescriptionIndex(localCellDescriptionsIndex), localCellDescriptionsIndex );
  auto& localCellDescriptions = getCellDescriptions(localCellDescriptionsIndex);

  int localElement = NotFound;
  for (unsigned int element = 0; element < localCellDescriptions.size(); ++element) {
    if ( localCellDescriptions[element].getSolverNumber()==receivedCellDescriptions[receivedElement].getSolverNumber() ) {
      localElement = element;
    }
    element++;
  }
  if ( localElement==NotFound ) { // We have already received data and allocated all memory in this case
    CellDescription& receivedCellDescription = receivedCellDescriptions[receivedElement];
    localCellDescriptions.push_back(receivedCellDescription);
    assertion1(receivedCellDescription.getType()==CellDescription::Type::Virtual,receivedCellDescription.toString());
    receivedCellDescriptions.erase(receivedCellDescriptions.begin()+receivedElement);
    return false;
  } else {
    CellDescription& receivedCellDescription = receivedCellDescriptions[receivedElement];
    CellDescription& localCellDescription = localCellDescriptions[localElement];
    localCellDescription.setParentType(receivedCellDescription.getType());

    if ( receivedCellDescriptions[receivedElement].getType()==CellDescription::Type::LeafProlongates ) {
      logDebug( "progressMeshRefinementInMergeWithWorker(...)","merging prolongated solution");

      assertion( localCellDescription.getType()==CellDescription::Type::Leaf ||
          localCellDescription.getType()==CellDescription::Type::Virtual);
      assertion(receivedCellDescription.getType()==CellDescription::Type::Leaf);
      assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getSolutionIndex()));
      assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getPreviousSolutionIndex()));

      // we know we have received data in this case
      localCellDescription.setType(CellDescription::Type::LeafProlongates);
      localCellDescription.setRefinementStatus(receivedCellDescription.getRefinementStatus());
      localCellDescription.setPreviousRefinementStatus(receivedCellDescription.getPreviousRefinementStatus());
      localCellDescription.setCommunicationStatus(LeafCommunicationStatus);
      localCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion

      ensureNecessaryMemoryIsAllocated(localCellDescription); // copy indices
      localCellDescription.setSolution(receivedCellDescription.getSolution());
      localCellDescription.setPreviousSolution(receivedCellDescription.getPreviousSolution());
      receivedCellDescription.setSolutionIndex(-1);
      receivedCellDescription.setPreviousSolutionIndex(-1);
      receivedCellDescription.setSolution(nullptr);
      receivedCellDescription.setPreviousSolution(nullptr);
      return true;
    } else {
      return false;
    }
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // send out data
  if ( cellDescription.getParentType()==CellDescription::Type::ParentCoarsens ) {
    sendDataToWorkerOrMasterDueToForkOrJoin(masterRank,cellDescriptionsIndex,element,
        peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy
  }
  eraseCellDescriptionIfNecessary(
      cellDescription,
      cellDescriptionsIndex, element,
      cellDescription.getParentType(), 0/*Not important here*/);
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  CellDescription& cellDescription = getCellDescription(localCellDescriptionsIndex,localElement);
  ensureFineGridCoarseGridConsistency(cellDescription,coarseGridCellDescriptionsIndex);
  cellDescription.setCreation(CellDescription::Creation::ReceivedFromWorker);

  CellDescription::Type parentType = CellDescription::Type::Erased;
  int parentAugmentationStatus     = 0;
  const int coarseGridElement = tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber());
  if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription = getCellDescription(cellDescription.getParentIndex(),coarseGridElement);

      if ( coarseGridCellDescription.getType()==CellDescription::Type::ParentCoarsens ) {
        mergeWithWorkerOrMasterDataDueToForkOrJoin(worker,localCellDescriptionsIndex,localElement,
            peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy
      }
      restrictVolumeDataIfParentCoarsens(cellDescription,coarseGridCellDescription);

      parentType               = coarseGridCellDescription.getType();
      parentAugmentationStatus = coarseGridCellDescription.getAugmentationStatus();
  }
  eraseCellDescriptionIfNecessary(cellDescription,localCellDescriptionsIndex,localElement,
      parentType,parentAugmentationStatus);

  // veto coarsening attempts of coarse grid cell
  if ( coarseGridElement != exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        cellDescription.getParentIndex(),coarseGridElement);
    vetoParentCoarseningRequestIfNecessary(cellDescription,coarseGridCellDescription);

    bool solverRequiresVerticalCommunication =
        ( isLeaf(cellDescription) || isVirtual(cellDescription) ) &&
        cellDescription.getAugmentationStatus()>MinimumAugmentationStatusForVirtualRefining;

    if (
      (coarseGridCellDescription.getType()==CellDescription::Type::Parent                   ||
      coarseGridCellDescription.getType()==CellDescription::Type::ParentRequestsCoarseningA ||
      coarseGridCellDescription.getType()==CellDescription::Type::ParentRequestsCoarseningB)
      &&
      solverRequiresVerticalCommunication
    ) {
       coarseGridCellDescription.setType(CellDescription::Type::ParentChecked);
    }
  }
  progressCoarseningOperationsInLeaveCell(cellDescription);

  decideOnRefinement(cellDescription);
}

///////////////////////////////////
// FORK OR JOIN
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
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

  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  if ( isLeaf(cellDescription) ) {
    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)",""
            "solution of solver " << cellDescription.getSolverNumber() << " sent to rank "<<toRank<<
                 ", cell: "<< x << ", level: " << level);

    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),
        cellDescriptionsIndex,cellDescription.toString());
    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()),
            cellDescriptionsIndex,cellDescription.toString());
    DataHeap::getInstance().sendData(
        static_cast<double*>(cellDescription.getSolution()),
        getDataPerCell(), toRank, x, level,messageType);
    DataHeap::getInstance().sendData(
        static_cast<double*>(cellDescription.getPreviousSolution()),
        getDataPerCell(), toRank, x, level,messageType);
  }
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  // allocate memory
  if ( messageType==peano::heap::MessageType::ForkOrJoinCommunication ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  } else if ( isLeaf(cellDescription) ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  }

  // receive data
  if ( isLeaf(cellDescription) ) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< x << ", level: " << level);

    getDataHeapEntries(cellDescription.getSolutionIndex()).clear();
    getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).clear();
    DataHeap::getInstance().receiveData(cellDescription.getSolutionIndex(),
        fromRank,x,level,messageType);
    DataHeap::getInstance().receiveData(cellDescription.getPreviousSolutionIndex(),
        fromRank,x,level,messageType);
  }
}

#endif
