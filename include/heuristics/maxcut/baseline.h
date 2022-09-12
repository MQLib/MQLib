#ifndef HEURISTICS_MAXCUT_BASELINE_H_
#define HEURISTICS_MAXCUT_BASELINE_H_

#include <vector>
#include "problem/max_cut_heuristic.h"

class Baseline : public MaxCutHeuristic {
 public:
  Baseline(const MaxCutInstance& mi, double runtime_limit, bool validation,
	   MaxCutCallback* mc);
};


#endif
