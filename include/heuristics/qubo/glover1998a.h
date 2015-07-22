#ifndef HEURISTICS_QUBO_GLOVER_1998A_H_
#define HEURISTICS_QUBO_GLOVER_1998A_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Glover1998aTabu {
 public:
  // Initialize to an empty tabu list for graph qi and recency list length t
  Glover1998aTabu(const QUBOInstance& qi);

  // Update the tabu lists to handle a new critical event solution -- remove the
  // oldest element in our recency list and add this one instead.
  void CriticalEvent(const QUBOSolution& x);

  // Get the tabu list values
  const std::vector<int>& get_tabuR() const {  return tabuR_;  }
  const std::vector<int>& get_tabuF() const {  return tabuF_;  }

 private:
  int t_;  // Length of adjacency circular list
  int N_;  // Number of variables in problem
  std::vector<int> tabuR_;  // Sum of t_ most recent critical solutions
  std::vector<int> tabuF_;  // Sum of all critical solutions
  std::vector<int> recent_;  // t_ most recent solutions as consec N_ elements
  int recent_pos_;  // The oldest position in recent_ (next to be removed)
};

class Glover1998aSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Glover1998aSolution
  Glover1998aSolution(const QUBOSolution &x);

  // Do the beginning parts of a phase (until critical event)
  void doPhaseBegin(int matchVal, int k, const Glover1998aTabu& tabu,
                    double PEN_R, double PEN_F);

  // Do the ending part of a phase (after critical event)
  void doPhaseEnd(int matchVal, int span);
};

class Glover1998a : public QUBOHeuristic {
 public:
  Glover1998a(const QUBOInstance& qi, double runtime_limit, bool validation,
              QUBOCallback *qc);

  // Update span variables to after one phase of an iteration
  void transferPhase();

 private:
  Glover1998aTabu tabu_;
  int k_;  // Number of iterations to use tabu list penalty
  int iter_span_;  // Number of iterations at current span value
  int span_;  // Number of moves after critical solution
  bool span_dir_increase_;  // Is the span increasing right now?
  int num_span_cycle_;  // How many span cycles have been completed?
  int N_;  // Number of variables
};

#endif
