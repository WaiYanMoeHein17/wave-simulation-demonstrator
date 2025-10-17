#include "ADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::ADERDGSolver::UpdateJob::UpdateJob(
  ADERDGSolver&                                      solver,
  CellDescription&                                   cellDescription,
  CellInfo&                                          cellInfo,
  const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers):
  tarch::multicore::jobs::Job(
      tarch::multicore::jobs::JobType::BackgroundTask,0,
      getHighPriorityTaskPriority()
  ), // ! always high priority
  _solver(solver),
  _cellDescription(cellDescription),
  _cellInfo(cellInfo),
  _boundaryMarkers(boundaryMarkers) {
  NumberOfReductionJobs.fetch_add(1);
}

bool exahype::solvers::ADERDGSolver::UpdateJob::run(bool runOnMasterThread) {
  _solver.updateBody(_cellDescription,_cellInfo,_boundaryMarkers);

  NumberOfReductionJobs.fetch_sub(1);
  assertion( NumberOfReductionJobs.load()>=0 );
  return false;
}


//
// @see UpdateJob
//
void exahype::solvers::ADERDGSolver::UpdateJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh    = static_cast<double*>(_cellDescription.getSolution());
  double* luhOld = static_cast<double*>(_cellDescription.getPreviousSolution());
  double* lduh   = static_cast<double*>(_cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(_cellDescription.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(_cellDescription.getFluctuation());

  _mm_prefetch(luh,    _MM_HINT_NTA);
  _mm_prefetch(luhOld, _MM_HINT_NTA);
  _mm_prefetch(lduh,   _MM_HINT_NTA);
  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}
