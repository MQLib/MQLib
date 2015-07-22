#ifndef HEURISTICS_MAXCUT_DESOUSA_2013_H_
#define HEURISTICS_MAXCUT_DESOUSA_2013_H_

#include "problem/max_cut_heuristic.h"

class deSousa2013 : public MaxCutHeuristic {
 public:
  deSousa2013(const MaxCutInstance &mi, double runtime_limit, bool validation,
              MaxCutCallback *mc);

 private:
  // Run through the algorithm, stopping when the p vector has converged.
  void AlgorithmIteration();
};

#endif
