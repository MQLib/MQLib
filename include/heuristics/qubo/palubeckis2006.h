#ifndef HEURISTICS_QUBO_PALUBECKIS_2006_H_
#define HEURISTICS_QUBO_PALUBECKIS_2006_H_

#include <vector>
#include "heuristics/qubo/palubeckis2004b.h"
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

// NOTE: In the paper Palubeckis2006, for a given solution x the authors
// describe a procedure in which the problem instance is adjusted (via a simple
// variable substitution) to one for which that solution is all 0s. This enables
// them to ignore the interaction terms between variables and only focus on
// the (adjusted) linear terms. After this adjustment, the linear terms are
// exactly equal to the diff_weights_ values maintained by the solution class.
// Therefore, instead of adjusting the problem instances, we just maintain the
// diff_weights_ as usual.

// Because a number of Palubeckis2006's methods are exactly the same as
// Palubeckis2004b's methods, we will extend Palubeckis2004bSolution and
// leverage it.
class Palubeckis2006Solution : public Palubeckis2004bSolution {
 public:
  // Convert from QUBOSolution to Palubeckis2006Solution
  Palubeckis2006Solution(const QUBOSolution &x);

  // Procedure GSP, which randomly perturbs a given solution with a preference
  // for high-quality moves. rTilde and b are parameters.
  void GSP(int rTilde, int b);

  // TS is the same as STS from the Palubeckis2004bSolution. best_objective
  // is a pointer to the best objective ever found in the master problem, and
  // mTilde is an iteration limit.
  void TS(double *best_objective, int mTilde) {  STS(best_objective, mTilde); }

  // LS is the same as LocalSearch from the Palubackis2004bSolution. z is a
  // pointer to the diff_weight_ check counter. LS does AllFirst1Swap updates.
  void LS(int *z)  {  LocalSearch(z); }
};

class Palubeckis2006 : public QUBOHeuristic {
 public:
  // ITS method from Palubeckis2006
  Palubeckis2006(const QUBOInstance& qi, double runtime_limit, bool validation,
                 QUBOCallback *qc);
};

#endif
