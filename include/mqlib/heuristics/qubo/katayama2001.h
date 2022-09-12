#ifndef HEURISTICS_QUBO_KATAYAMA_2001_H_
#define HEURISTICS_QUBO_KATAYAMA_2001_H_

#include "mqlib/heuristics/qubo/qubo_solution.h"
#include "mqlib/problem/qubo_heuristic.h"
#include "mqlib/problem/qubo_instance.h"

namespace mqlib {

    class Katayama2001Solution : public QUBOSolution {
    public:
        // Convert from QUBOSolution to Katayama2001
        Katayama2001Solution(const QUBOSolution &x);

        // Like the swap in local search, but has a chance of happening
        // even if worse
        bool SASwap(int k, double T);
    };

    class Katayama2001 : public QUBOHeuristic {
    public:
        Katayama2001(const QUBOInstance &qi, double runtime_limit, bool validation,
                     QUBOCallback *qc);
    };


}

#endif
