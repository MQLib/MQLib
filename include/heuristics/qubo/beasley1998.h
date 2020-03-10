#ifndef HEURISTICS_QUBO_BEASLEY_1998_H_
#define HEURISTICS_QUBO_BEASLEY_1998_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Beasley1998Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Beasley1998
  Beasley1998Solution(const QUBOSolution &x);

  void SA(double T);
  void LocalSearch(int &t);
  int TS(std::vector<int> &L, int iter, double vStar, int &t);
};

// Simulated Annealing algorithm
class Beasley1998SA : public QUBOHeuristic {
 public:
  Beasley1998SA(const QUBOInstance& qi, double runtime_limit, bool validation,
                QUBOCallback *qc);
};

// Tabu Search algorithm
class Beasley1998TS : public QUBOHeuristic {
 public:
  Beasley1998TS(const QUBOInstance& qi, double runtime_limit, bool validation,
                QUBOCallback *qc);
};

#endif
