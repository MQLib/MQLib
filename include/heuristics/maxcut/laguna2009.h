#ifndef HEURISTICS_MAXCUT_LAGUNA_2009_H_
#define HEURISTICS_MAXCUT_LAGUNA_2009_H_

#include <vector>
#include "heuristics/maxcut/max_cut_solution.h"
#include "problem/max_cut_heuristic.h"

class Laguna2009CE : public MaxCutHeuristic {
 public:
  Laguna2009CE(const MaxCutInstance& mi, double runtime_limit, bool validation,
               MaxCutCallback *mc);
};

class Laguna2009HCE : public MaxCutHeuristic {
 public:
  Laguna2009HCE(const MaxCutInstance& mi, double runtime_limit, bool validation,
                MaxCutCallback *mc);

 private:
  // Run the algorithm in Fig. 2 on the population X
  void LocallyOptimize(std::vector<FirstFixedMaxCutSolution>* X);
};

#endif
