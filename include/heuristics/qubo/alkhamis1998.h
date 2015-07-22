#ifndef HEURISTICS_QUBO_ALKHAMIS_1998_H_
#define HEURISTICS_QUBO_ALKHAMIS_1998_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

// alkhamis1998 proposed a simulated annealing approach very similar to the one
// in the appendix of korst1989; it differs mainly in the parameters selected
// for the SA.

class Alkhamis1998Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Alkhamis1998Solution
 Alkhamis1998Solution(const QUBOSolution &x) :
  QUBOSolution(x) {}

  // Run the full SA algorithm from section 2.2 on this solution.
  void SA(double T_initial, int iteration);
};

// Repeated sumulation annealing from random initial solutions
class Alkhamis1998 : public QUBOHeuristic {
 public:
  Alkhamis1998(const QUBOInstance& qi, double runtime_limit, bool validation,
               QUBOCallback *qc);
};

#endif
