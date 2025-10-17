#include "exahype/solvers/LimitingADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  LimitingADERDGSolver&                              solver,
  SolverPatch&                                       solverPatch,
  CellInfo&                                          cellInfo,
  const double                                       predictionTimeStamp,
  const double                                       predictionTimeStepSize,
  const bool                                         isFirstTimeStepOfBatch,
  const bool                                         isLastTimeStepOfBatch,
  const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
  const bool                                         isSkeletonJob):
  tarch::multicore::jobs::Job(
      tarch::multicore::jobs::JobType::BackgroundTask,0,
      getTaskPriority(isLastTimeStepOfBatch)
  ),
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _predictionTimeStamp   (predictionTimeStamp),
  _predictionTimeStepSize(predictionTimeStepSize),
  _isFirstTimeStepOfBatch(isFirstTimeStepOfBatch),
  _isLastTimeStepOfBatch (isLastTimeStepOfBatch),
  _boundaryMarkers(boundaryMarkers),
  _isSkeletonJob(isSkeletonJob) {
  NumberOfReductionJobs.fetch_add(1);
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    NumberOfEnclaveJobs.fetch_add(1);
  }
}

bool exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::run(bool runOnMasterThread) {
  _solver.fusedTimeStepBody(
      _solverPatch,_cellInfo,
      _predictionTimeStamp,_predictionTimeStepSize,
      _isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
      _boundaryMarkers,_isSkeletonJob,false/*mustBeDoneImmedetially*/);

  NumberOfReductionJobs.fetch_sub(1);
  assertion( NumberOfReductionJobs.load()>=0 );
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_sub(1);
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs.fetch_sub(1);
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
  return false;
}
