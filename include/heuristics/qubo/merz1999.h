#ifndef HEURISTICS_QUBO_MERZ_1999_H_
#define HEURISTICS_QUBO_MERZ_1999_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Merz1999Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Merz1999Solution
  Merz1999Solution(const QUBOSolution &x);

  // Crossover (HUX) and Local Search combined
  static Merz1999Solution HUX(const QUBOInstance& qi,
                              const Merz1999Solution& parent_a,
                              const Merz1999Solution& parent_b,
                              QUBOHeuristic *heuristic) {
    return Merz1999Solution(qi, parent_a, parent_b, heuristic);
  }

  // Standard uniform crossover
  static Merz1999Solution Crossover(const QUBOInstance& qi,
                                    const Merz1999Solution& parent_a,
                                    const Merz1999Solution& parent_b,
                                    QUBOHeuristic *heuristic) {
    return Merz1999Solution(qi, parent_a, parent_b, heuristic, 0.5);
  }

  void RestartMutate();
  void Mutate();

 private:
  Merz1999Solution(const QUBOInstance& qi,
                   const Merz1999Solution& parent_a,
                   const Merz1999Solution& parent_b,
                   QUBOHeuristic *heuristic);

  Merz1999Solution(const QUBOInstance& qi,
                   const Merz1999Solution& parent_a,
                   const Merz1999Solution& parent_b,
                   QUBOHeuristic *heuristic,
                   double dummy);
};


class Merz1999 : public QUBOHeuristic {
 public:
  Merz1999(const QUBOInstance& qi, double runtime_limit, bool validation,
	   QUBOCallback *qc, int version);
};

// The local search variant
class Merz1999GLS : public Merz1999 {
 public:
 Merz1999GLS(const QUBOInstance& qi, double runtime_limit, bool validation,
             QUBOCallback *qc) :
  Merz1999(qi, runtime_limit, validation, qc, 0) {}
};

// The crossover-only variant
class Merz1999Crossover : public Merz1999 {
 public:
 Merz1999Crossover(const QUBOInstance& qi, double runtime_limit,
		   bool validation, QUBOCallback *qc) :
  Merz1999(qi, runtime_limit, validation, qc, 1) {}
};

// The mutation-only variant
class Merz1999Mutation : public Merz1999 {
 public:
 Merz1999Mutation(const QUBOInstance& qi, double runtime_limit,
		  bool validation, QUBOCallback *qc) :
  Merz1999(qi, runtime_limit, validation, qc, 2) {}
};

#endif
