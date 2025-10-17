#include "LimitingADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::LimitingADERDGSolver::UpdateJob::UpdateJob(
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

bool exahype::solvers::LimitingADERDGSolver::UpdateJob::run(bool runOnMasterThread) {
  _solver.updateBody(_solverPatch,_cellInfo,_boundaryMarkers);

  NumberOfReductionJobs.fetch_sub(1);
  assertion( NumberOfReductionJobs.load()>=0 );

  return false;
}

//
// @see UpdateJob
//
void exahype::solvers::LimitingADERDGSolver::UpdateJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh    = static_cast<double*>(_solverPatch.getSolution());
  double* luhOld = static_cast<double*>(_solverPatch.getPreviousSolution());
  double* lduh   = static_cast<double*>(_solverPatch.getUpdate());
  double* lQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(_solverPatch.getFluctuation());

  _mm_prefetch(luh,    _MM_HINT_NTA);
  _mm_prefetch(luhOld, _MM_HINT_NTA);
  _mm_prefetch(lduh,   _MM_HINT_NTA);
  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}
