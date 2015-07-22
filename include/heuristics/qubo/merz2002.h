#ifndef HEURISTICS_QUBO_MERZ_2002_H_
#define HEURISTICS_QUBO_MERZ_2002_H_

#include "heuristics/qubo/qubo_partial_solution.h"
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Merz2002PartialSolution : public QUBOPartialSolution {
 public:
  // Create a "greedy random" solution using procedure RandomizedGreedy
  static Merz2002PartialSolution RandomizedGreedy(const QUBOInstance& qi,
						  QUBOHeuristic *heuristic) {
    return Merz2002PartialSolution(qi, heuristic);
  }

 private:
  // Create a "greedy random" solution using procedure RandomizedGreedy
  Merz2002PartialSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic);
};

class Merz2002QUBOSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Merz2002Solution
 Merz2002QUBOSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

  // Perform k-opt using procedure Local-Search-k-opt
  void KOpt();
};

// Base class: specified by random vs. greedy construction as well as
// 1-opt vs. k-opt vs. no local search
class Merz2002 : public QUBOHeuristic {
 public:
  enum LStype {NONE, ONEOPT, KOPT};

  Merz2002(const QUBOInstance& qi, double runtime_limit, bool validation,
	   QUBOCallback *qc, bool greedy, LStype type);
};

// The greedy construction and no local search variant
class Merz2002Greedy : public Merz2002 {
 public:
 Merz2002Greedy(const QUBOInstance& qi, double runtime_limit, bool validation,
                QUBOCallback *qc) :
  Merz2002(qi, runtime_limit, validation, qc, true, NONE) {}
};

// The random construction and 1-opt variant
class Merz2002OneOpt : public Merz2002 {
 public:
 Merz2002OneOpt(const QUBOInstance& qi, double runtime_limit, bool validation,
                QUBOCallback *qc) :
  Merz2002(qi, runtime_limit, validation, qc, false, ONEOPT) {}
};

// The random construction and k-opt variant
class Merz2002KOpt : public Merz2002 {
 public:
 Merz2002KOpt(const QUBOInstance& qi, double runtime_limit, bool validation,
              QUBOCallback *qc) :
  Merz2002(qi, runtime_limit, validation, qc, false, KOPT) {}
};

// The greedy construction and k-opt variant
class Merz2002GreedyKOpt : public Merz2002 {
 public:
 Merz2002GreedyKOpt(const QUBOInstance& qi, double runtime_limit,
		    bool validation, QUBOCallback *qc) :
  Merz2002(qi, runtime_limit, validation, qc, true, KOPT) {}
};

#endif
