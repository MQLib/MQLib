#ifndef HEURISTICS_MAXCUT_BASELINE_H_
#define HEURISTICS_MAXCUT_BASELINE_H_

#include <vector>
#include "mqlib/problem/max_cut_heuristic.h"

namespace mqlib {

    class Baseline : public MaxCutHeuristic {
    public:
        Baseline(const MaxCutInstance &mi, double runtime_limit, bool validation,
                 MaxCutCallback *mc);
    };
}

#endif
