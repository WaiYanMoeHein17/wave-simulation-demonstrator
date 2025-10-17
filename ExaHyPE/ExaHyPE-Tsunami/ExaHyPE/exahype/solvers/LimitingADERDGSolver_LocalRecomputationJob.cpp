#include "LimitingADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::LimitingADERDGSolver::LocalRecomputationJob::LocalRecomputationJob(
  LimitingADERDGSolver&                              solver,
  SolverPatch&                                       solverPatch,
  CellInfo&                                          cellInfo,
  const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers):
  tarch::multicore::jobs::Job(
      tarch::multicore::jobs::JobType::BackgroundTask,0,
      getHighPriorityTaskPriority()
  ), // ! always high priority
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _boundaryMarkers(boundaryMarkers) {
  NumberOfReductionJobs.fetch_add(1);
}

bool exahype::solvers::LimitingADERDGSolver::LocalRecomputationJob::run(bool runOnMasterThread) {
  double admissibleTimeStepSize =
      _solver.localRecomputationBody(_solverPatch,_cellInfo,_boundaryMarkers);

  _solver.updateAdmissibleTimeStepSize(admissibleTimeStepSize);

  NumberOfReductionJobs.fetch_sub(1);
  assertion( NumberOfReductionJobs.load()>=0 );

  return false;
}


//
// @see UpdateJob
//
void exahype::solvers::LimitingADERDGSolver::LocalRecomputationJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  if (
      _solver.isInvolvedInLocalRecomputation(_solverPatch) &&
      _solverPatch.getRefinementStatus()>=_solver._solver->_minRefinementStatusForTroubledCell-1
  ) {
    LimiterPatch& limiterPatch = _solver.getLimiterPatch(_solverPatch,_cellInfo);
    double* luh    = static_cast<double*>(limiterPatch.getSolution());
    double* luhOld = static_cast<double*>(limiterPatch.getPreviousSolution());

    _mm_prefetch(luh,    _MM_HINT_NTA);
    _mm_prefetch(luhOld, _MM_HINT_NTA);
  }
  #endif
}
