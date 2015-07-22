#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/maxcut/max_cut_simple_solution.h"
#include "util/random.h"

// Empty solution
MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutInstance& mi,
                                           MaxCutHeuristic *heuristic,
                                           int initVal) :
  BaseSolution(mi.get_size(), initVal),
  mi_(mi),
  heuristic_(heuristic) {}

// Solution with provided assignments and weight
MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutInstance& mi,
                                           MaxCutHeuristic *heuristic,
                                           const std::vector<int>& assignments,
                                           double weight) :
  BaseSolution(assignments, weight),
  mi_(mi),
  heuristic_(heuristic) {}

// Random solution (p=0.5 for each vertex)
MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutInstance& mi,
                                           MaxCutHeuristic *heuristic,
                                           int ignored1, int ignored2):
  BaseSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // Random assignments in {-1, 1}
  for (int i=0; i < mi.get_size(); ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_
  PopulateFromAssignments();
}

// Random solution (p provided for each vertex)
MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutInstance& mi,
                                           const std::vector<double>& p,
                                           MaxCutHeuristic *heuristic) :
  BaseSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // Random assignments in {-1, 1} using probabilities p of being set 1
  for (int i=0; i < mi.get_size(); ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_
  PopulateFromAssignments();
}

MaxCutSimpleSolution::MaxCutSimpleSolution(const QUBOSimpleSolution& sol,
                                           const MaxCutInstance& mi,
                                           MaxCutHeuristic *heuristic) :
  BaseSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // To map from QUBO to Max-Cut, the size stays the same and the assignments
  // map 0 -> -1 and 1 -> 1. The objective stays the same. We have initialized
  // to a solution with all -1 values and objective 0, so we need to
  // change to assignment 1 when appropriate and copy the objective value.
  const std::vector<int>& qubo_assignments = sol.get_assignments();
  if (qubo_assignments.size() != N_) {
    std::cout << "ERROR: Instance size mismatch when building " <<
      "MaxCutSimpleSolution from QUBOSimpleSolution" << std::endl;
    exit(1);
  }
  for (int i=0; i < N_; ++i) {
    if (qubo_assignments[i] == 1) {
      assignments_[i] = 1;
    }
  }
  weight_ = sol.get_weight();
}

MaxCutSimpleSolution& MaxCutSimpleSolution::operator=(const MaxCutSimpleSolution &rhs) {
  BaseSolution::operator=(rhs);
  // Can't copy over mi_ because it's a reference
  heuristic_ = rhs.heuristic_;
  return *this;
}

MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutSimpleSolution &x) :
  BaseSolution(x),
  mi_(x.mi_),
  heuristic_(x.heuristic_) {}

void MaxCutSimpleSolution::PopulateFromAssignments() {
  weight_ = 0.0;

  for (auto it = mi_.get_all_edges_begin(); it != mi_.get_all_edges_end();
       ++it) {
    if (assignments_[it->first.first] != assignments_[it->first.second]) {
      // These two nodes are in different sets
      weight_ += it->second;
    }
  }
}
