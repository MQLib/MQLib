#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_instance.h"
#include "util/random.h"

// Empty solution
QUBOSolution::QUBOSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic) :
  ExtendedSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {}

// All-zero solution
QUBOSolution::QUBOSolution(const QUBOInstance &qi, QUBOHeuristic *heuristic,
                           int ignored1, int ignored2) :
  ExtendedSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {
  diff_weights_ = qi.get_lin();
}

// Random solution (p=0.5 for each vertex)
QUBOSolution::QUBOSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic,
			   int ignored):
  ExtendedSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {
  // Random assignments in {0, 1}
  for (int i=0; i < qi.get_size(); ++i) {
    assignments_[i] = Random::RandInt(0, 1);
  }

  // Obtain weight_ and diff_weights_
  PopulateFromAssignments();
}

// Random solution (p provided for each variable)
QUBOSolution::QUBOSolution(const QUBOInstance& qi,
			   const std::vector<double>& p,
			   QUBOHeuristic *heuristic) :
  ExtendedSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {
  // Random assignments in {0, 1} using probabilities p of being set 1
  for (int i=0; i < qi.get_size(); ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : 0;
  }

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Initialize from vector of assignments
QUBOSolution::QUBOSolution(const std::vector<int>& assignments,
                           const QUBOInstance& qi,
                           QUBOHeuristic *heuristic) :
  ExtendedSolution(qi.get_size(), 0),
  qi_(qi),
  heuristic_(heuristic) {
  // Copy over assignments
  assignments_ = assignments;

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Does all the updates required when the node at update_index has its cut set
// inclusion flipped. x is updated to show this flipped inclusion, and the
// difference in weights to own set and other set are updated.
void QUBOSolution::UpdateCutValues(int update_index, std::vector<int>* x,
				   std::vector<double>* diff_weights,
				   double *objective) const {
  *objective += (*diff_weights)[update_index];
  (*x)[update_index] = 1 - (*x)[update_index];
  (*diff_weights)[update_index] = -(*diff_weights)[update_index];

  // Iterate the set of all neighbors for node update_index
  for (auto iter = qi_.get_nonzero_begin(update_index);
       iter != qi_.get_nonzero_end(update_index); ++iter) {
    // var update_index has non-zero interaction score q_ij with var j
    int j = iter->first;
    double q_ij = iter->second;
    if ((*x)[update_index] == (*x)[j]) {
      (*diff_weights)[j] -= 2.0 * q_ij;
    } else {
      (*diff_weights)[j] += 2.0 * q_ij;
    }
  }
}

void QUBOSolution::AllBest2Swap() {
  // Take all profitable 2-moves, taking the most profitable first.
  while (1) {
    double best_move = 0.0;
    int best_i = -1;
    int best_j = -1;
    for (auto iter = qi_.get_all_nonzero_begin();
	 iter != qi_.get_all_nonzero_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double q_ij = iter->second;
      // Benefit DW_i + DW_j + 2q_ij if started with same assignment
      // Benefit DW_i + DW_j - 2q_ij if started with different assignment
      double benefit = diff_weights_[i] + diff_weights_[j] +
	(4 * (assignments_[i] == assignments_[j]) - 2) * q_ij;
      if (benefit > best_move) {
	best_move = benefit;
	best_i = i;
	best_j = j;
      }
    }
    if (best_i < 0 || !ImprovingMove(best_move)) {
      // No more profitable moves
      break;
    }

    // Update the diff_weights_ variables and objective
    UpdateCutValues(best_i);
    UpdateCutValues(best_j);
  }
}

void QUBOSolution::AllFirst2Swap() {
  // Take all profitable 2-moves, taking the first one we find when scanning
  // the non-zero q_ij values sequentially.
  int move_made = 1;
  while (move_made) {
    move_made = 0;
    for (auto iter = qi_.get_all_nonzero_begin();
	 iter != qi_.get_all_nonzero_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double q_ij = iter->second;
      // Benefit DW_i + DW_j + 2q_ij if started with same assignment
      // Benefit DW_i + DW_j - 2q_ij if started with different assignment
      double benefit = diff_weights_[i] + diff_weights_[j] +
	(4 * (assignments_[i] == assignments_[j]) - 2) * q_ij;
      if (ImprovingMove(benefit)) {
	UpdateCutValues(i);
	UpdateCutValues(j);
	move_made = 1;
	break;
      }
    }
  }
}

QUBOSolution& QUBOSolution::operator=(const QUBOSolution &rhs) {
  ExtendedSolution::operator=(rhs);
  heuristic_ = rhs.heuristic_;
  return *this;
}

QUBOSolution::QUBOSolution(const QUBOSolution &x) :
  ExtendedSolution(x),
  qi_(x.qi_),
  heuristic_(x.heuristic_) {}

void QUBOSolution::PopulateFromAssignments() {
  weight_ = 0.0;
  diff_weights_.assign(N_, 0.0);

  // First, deal with the linear terms based on the assignments
  for (int i=0; i < N_; ++i) {
    if (assignments_[i]) {
      diff_weights_[i] -= qi_.get_lin()[i];
      weight_ += qi_.get_lin()[i];
    } else {
      diff_weights_[i] += qi_.get_lin()[i];
    }
  }

  // Next, deal with the pairs of variables with non-zero interaction terms
  for (auto it = qi_.get_all_nonzero_begin(); it != qi_.get_all_nonzero_end();
       ++it) {
    // There is a non-zero interaction between variables i and j, with weight
    // q_ij.
    int i = it->first.first;
    int j = it->first.second;
    double q_ij = it->second;

    if (assignments_[i] && assignments_[j]) {
      // Both variables are set to 1
      weight_ += 2.0 * q_ij;
      diff_weights_[i] -= 2.0 * q_ij;
      diff_weights_[j] -= 2.0 * q_ij;
    } else if (assignments_[i] && !assignments_[j]) {
      // Only i is set
      diff_weights_[j] += 2.0 * q_ij;
    } else if (!assignments_[i] && assignments_[j]) {
      // Only j is set
      diff_weights_[i] += 2.0 * q_ij;
    }
  }
}

// Copy over a QUBOPartialSolution to be a full QUBOSolution
QUBOSolution::QUBOSolution(const QUBOPartialSolution& x) :
  ExtendedSolution(x.get_qi().get_size(), 0),
  qi_(x.get_qi()),
  heuristic_(x.get_heuristic()) {
  if (x.get_num_frac() > 0) {
    std::cout << "Cannot copy over fractional QUBOPartialSolution" << std::endl;
    exit(0);
  }

  // Copy over objective from the partial solution
  weight_ = x.get_weight();

  // Copy over the assignments from the partial solution, converting to int.
  // Copy over the appropriate diff_weights_ value, as well.
  for (int i=0; i < N_; ++i) {
    assignments_[i] = (int)(x.get_assignments()[i]);
    if (assignments_[i]) {
      diff_weights_[i] = x.get_diff0()[i];
    } else {
      diff_weights_[i] = x.get_diff1()[i];
    }
  }
}
