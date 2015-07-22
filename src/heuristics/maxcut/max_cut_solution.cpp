#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/maxcut/max_cut_solution.h"
#include "util/random.h"

// Empty solution
MaxCutSolution::MaxCutSolution(const MaxCutInstance& mi,
			       MaxCutHeuristic *heuristic,
			       int initVal) :
  ExtendedSolution(mi.get_size(), initVal),
  mi_(mi),
  heuristic_(heuristic) {}

// Random solution (p=0.5 for each vertex)
MaxCutSolution::MaxCutSolution(const MaxCutInstance& mi,
			       MaxCutHeuristic *heuristic, int ignored1,
			       int ignored2):
  ExtendedSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // Random assignments in {-1, 1}
  for (int i=0; i < mi.get_size(); ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_ and diff_weights_
  PopulateFromAssignments();
}

// Random solution (p provided for each vertex)
MaxCutSolution::MaxCutSolution(const MaxCutInstance& mi,
			       const std::vector<double>& p,
			       MaxCutHeuristic *heuristic) :
  ExtendedSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // Random assignments in {-1, 1} using probabilities p of being set 1
  for (int i=0; i < mi.get_size(); ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Initialize from vector of assignments
MaxCutSolution::MaxCutSolution(const std::vector<int>& assignments,
                               const MaxCutInstance& mi,
                               MaxCutHeuristic *heuristic) :
  ExtendedSolution(mi.get_size(), -1),
  mi_(mi),
  heuristic_(heuristic) {
  // Copy over assignments
  assignments_ = assignments;

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Does all the updates required when the node at update_index has its cut set
// inclusion flipped. x is updated to show this flipped inclusion, and the
// difference in weights to own set and other set are updated.
void MaxCutSolution::UpdateCutValues(int update_index, std::vector<int>* x,
				     std::vector<double>* diff_weights,
				     double *objective) const {
  *objective += (*diff_weights)[update_index];
  (*x)[update_index] = -(*x)[update_index];
  (*diff_weights)[update_index] = -(*diff_weights)[update_index];

  // Iterate the set of all neighbors for node update_index
  for (auto iter = mi_.get_edges_begin(update_index);
       iter != mi_.get_edges_end(update_index); ++iter) {
    int j = iter->first;  // There is an edge (update_index, j) in the graph
    double w_ij = iter->second;  // The edge has weight w_ij
    (*diff_weights)[j] += 2.0 * (*x)[update_index] * (*x)[j] * w_ij;
  }
}

void MaxCutSolution::AllBest2Swap(int startpos) {
  // Take all profitable 2-moves, taking the most profitable first.
  while (true) {
    double best_move = 0.0;
    int best_i = -1;
    int best_j = -1;
    for (auto iter = mi_.get_all_edges_begin();
	 iter != mi_.get_all_edges_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double w_ij = iter->second;
      double benefit = diff_weights_[i] + diff_weights_[j] -
	2.0 * assignments_[i] * assignments_[j] * w_ij;
      // Only note move if it improves on best found so far and i and j
      // are both large enough. Note that there's limited efficiency penalty
      // compared to writing two separate functions (one without startpos and
      // one with it) because the last two conditions are only checked
      // the few times when we've found an improving 2-move. The same argument
      // also applies for AllFirst2Swap() below.
      if (benefit > best_move && i >= startpos && j >= startpos) {
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

void MaxCutSolution::AllFirst2Swap(int startpos) {
  // Take all profitable 2-moves, taking the first one we find when scanning
  // the edges sequentially.
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (auto iter = mi_.get_all_edges_begin();
	 iter != mi_.get_all_edges_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double w_ij = iter->second;
      double benefit = diff_weights_[i] + diff_weights_[j] -
	2.0 * assignments_[i] * assignments_[j] * w_ij;
      if (ImprovingMove(benefit) && i >= startpos && j >= startpos) {
	UpdateCutValues(i);
	UpdateCutValues(j);
	move_made = true;
	break;
      }
    }
  }
}

MaxCutSolution& MaxCutSolution::operator=(const MaxCutSolution &rhs) {
  ExtendedSolution::operator=(rhs);
  // Can't copy over mi_ because it's a reference
  heuristic_ = rhs.heuristic_;
  return *this;
}

MaxCutSolution::MaxCutSolution(const MaxCutSolution &x) :
  ExtendedSolution(x),
  mi_(x.mi_),
  heuristic_(x.heuristic_) {}

void MaxCutSolution::PopulateFromAssignments() {
  weight_ = 0.0;
  diff_weights_.assign(N_, 0.0);

  for (auto it = mi_.get_all_edges_begin(); it != mi_.get_all_edges_end();
       ++it) {
    if (assignments_[it->first.first] == assignments_[it->first.second]) {
      // These two nodes are in the same set
      diff_weights_[it->first.first] += it->second;
      diff_weights_[it->first.second] += it->second;
    } else {
      // These two nodes are in different sets
      weight_ += it->second;
      diff_weights_[it->first.first] -= it->second;
      diff_weights_[it->first.second] -= it->second;
    }
  }
}

// Empty solution
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(const MaxCutInstance& mi,
						   MaxCutHeuristic *heuristic,
						   int fixedVal) :
  MaxCutSolution(mi, heuristic, fixedVal),
  fixedVal_(fixedVal) {}

// Random solution
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(const MaxCutInstance& mi,
						   MaxCutHeuristic *heuristic,
						   int fixedVal, int ignored2):
  MaxCutSolution(mi, heuristic, fixedVal),
  fixedVal_(fixedVal) {
  // Random assignments in {-1, 1} (don't assign the first one; already handled
  // by the MaxCutSolution initializer).
  for (int i=1; i < mi.get_size(); ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_ and diff_weights_
  PopulateFromAssignments();
}

// Random solution from vertex probabilities
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(const MaxCutInstance& mi,
						   const std::vector<double>& p,
						   MaxCutHeuristic *heuristic,
						   int fixedVal) :
  MaxCutSolution(mi, heuristic, fixedVal),
  fixedVal_(fixedVal) {
  // Random assignments in {-1, 1} using probabilities p of being set 1 (don't
  // assign the first one; already handled by MaxCutSolution initializer)
  for (int i=1; i < mi.get_size(); ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Does all the updates required when the node at update_index has its cut set
// inclusion flipped. x is updated to show this flipped inclusion, and the
// difference in weights to own set and other set are updated.
void FirstFixedMaxCutSolution::UpdateCutValues(int update_index,
					       std::vector<int>* x,
					       std::vector<double>* diff_wts,
					       double *objective) const {
  if (update_index == 0) {
    std::cout << "Error: flipping first index of a FirstFixedMaxCutSolution" <<
      std::endl;
    exit(0);
  }
  MaxCutSolution::UpdateCutValues(update_index, x, diff_wts, objective);
}

FirstFixedMaxCutSolution&
FirstFixedMaxCutSolution::operator=(const FirstFixedMaxCutSolution &rhs) {
  MaxCutSolution::operator=(rhs);
  fixedVal_ = rhs.fixedVal_;
  return *this;
}

FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(const
						   FirstFixedMaxCutSolution &x)
  : MaxCutSolution(x),
    fixedVal_(x.fixedVal_) {}

void FirstFixedMaxCutSolution::PopulateFromAssignments() {
  if (assignments_[0] != fixedVal_) {
    std::cout << "Error: wrong start val in PopulateFromAssignments" <<
      std::endl;
    exit(0);
  }
  MaxCutSolution::PopulateFromAssignments();
}

// Copy over a MaxCutPartialSolution to be a full MaxCutSolution
MaxCutSolution::MaxCutSolution(const MaxCutPartialSolution& x) :
  ExtendedSolution(x.get_mi().get_size(), 0),
  mi_(x.get_mi()),
  heuristic_(x.get_heuristic()) {
  if (x.get_num_unassigned() > 0) {
    std::cout << "Cannot copy over MaxCutPartialSolution with unassigned nodes"
	      << std::endl;
    exit(0);
  }

  // Copy over assignments and objective value from partial solution
  assignments_.assign(x.get_assignments().begin(), x.get_assignments().end());
  weight_ = x.get_weight();

  // Assign diff_weights_ from gainS and gainNS of the partial solution
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] == 1) {
      diff_weights_[i] = x.get_gainNS()[i];
    } else {
      diff_weights_[i] = x.get_gainS()[i];
    }
  }
}
