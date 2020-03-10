#include "heuristics/maxcut/baseline.h"
#include "heuristics/maxcut/max_cut_solution.h"

Baseline::Baseline(const MaxCutInstance& mi, double runtime_limit,
		   bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  int iters = static_cast<int>(runtime_limit);
  for (int iter=0; ; ++iter) {
    MaxCutSolution x = MaxCutSolution::RandomSolution(mi, this);
    x.AllFirst1Swap();
    if (iter == iters - 1) {
      Report(x);
      return;
    }
  }
}
