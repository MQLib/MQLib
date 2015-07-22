#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/qubo/qubo_simple_solution.h"
#include "util/random.h"

// Empty solution
QUBOSimpleSolution::QUBOSimpleSolution(const QUBOInstance& qi,
                                       QUBOHeuristic *heuristic,
                                       int initVal) :
  BaseSolution(qi.get_size(), initVal),
  qi_(qi),
  heuristic_(heuristic) {}

// Solution with provided assignments and weights
QUBOSimpleSolution::QUBOSimpleSolution(const QUBOInstance& qi,
                                       QUBOHeuristic *heuristic,
                                       const std::vector<int>& assignments,
                                       double weight) :
  BaseSolution(assignments, weight),
  qi_(qi),
  heuristic_(heuristic) {}

// Random solution (p=0.5 for each vertex)
QUBOSimpleSolution::QUBOSimpleSolution(const QUBOInstance& qi,
                                       QUBOHeuristic *heuristic,
                                       int ignored1, int ignored2):
  BaseSolution(qi.get_size(), -1),
  qi_(qi),
  heuristic_(heuristic) {
  // Random assignments in {0, 1}
  for (int i=0; i < qi.get_size(); ++i) {
    assignments_[i] = Random::RandInt(0, 1);
  }

  // Obtain weight_
  PopulateFromAssignments();
}

// Random solution (p provided for each vertex)
QUBOSimpleSolution::QUBOSimpleSolution(const QUBOInstance& qi,
                                       const std::vector<double>& p,
                                       QUBOHeuristic *heuristic) :
  BaseSolution(qi.get_size(), -1),
  qi_(qi),
  heuristic_(heuristic) {
  // Random assignments in {0, 1} using probabilities p of being set 1
  for (int i=0; i < qi.get_size(); ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : 0;
  }

  // Obtain weights_
  PopulateFromAssignments();
}

QUBOSimpleSolution::QUBOSimpleSolution(const MaxCutSimpleSolution& sol,
                                       const QUBOInstance& qi,
                                       QUBOHeuristic *heuristic) :
  BaseSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {
  // To map from Max-Cut to QUBO, the size decreases by 1 because we are
  // removing the master node (the last node) that we added when we reduced
  // QUBO to Max-Cut. The assigned values for the variables are 1 if their
  // associated node was across the cut from the master node and 0 if their
  // associated node was on the same side of the cut as the master node.
  // The objective stays the same. We have initialized to a solution with all
  // 0 values and objective 0, so we need to change to assignment 1 when
  // appropriate and copy the objective value.
  const std::vector<int>& mc_assignments = sol.get_assignments();
  if (mc_assignments.size() != N_+1) {
    std::cout << "ERROR: Instance size mismatch when building " <<
      "QUBOSimpleSolution from MaxCutSimpleSolution" << std::endl;
    exit(1);
  }

  for (int i=0; i < N_; ++i) {
    if (mc_assignments[i] != mc_assignments[N_]) {
      assignments_[i] = 1;
    }
  }
  weight_ = sol.get_weight();
}

QUBOSimpleSolution& QUBOSimpleSolution::operator=(const QUBOSimpleSolution &rhs) {
  BaseSolution::operator=(rhs);
  // Can't copy over qi_ because it's a reference
  heuristic_ = rhs.heuristic_;
  return *this;
}

QUBOSimpleSolution::QUBOSimpleSolution(const QUBOSimpleSolution &x) :
  BaseSolution(x),
  qi_(x.qi_),
  heuristic_(x.heuristic_) {}

void QUBOSimpleSolution::PopulateFromAssignments() {
  weight_ = 0.0;

  // Linear terms
  const std::vector<double>& lin = qi_.get_lin();
  for (int i=0; i < N_; ++i) {
    weight_ += assignments_[i] * lin[i];
  }

  // Interaction terms
  for (auto it = qi_.get_all_nonzero_begin(); it != qi_.get_all_nonzero_end();
       ++it) {
    if (assignments_[it->first.first] == 1 &&
        assignments_[it->first.second] == 1) {
      // Both nodes are set, so increment the objective
      weight_ += 2.0 * it->second;
    }
  }
}
