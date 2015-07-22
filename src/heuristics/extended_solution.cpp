#include <algorithm>
#include <math.h>
#include "heuristics/extended_solution.h"

ExtendedSolution::ExtendedSolution(int N, int init_assignment) :
  BaseSolution(N, init_assignment),
  diff_weights_(N, 0.0) {}

void ExtendedSolution::AllBest1Swap(int startpos) {
  // Take all profitable 1-moves, taking the most profitable first.
  while (true) {
    double best_move = 0.0;
    int best_pos = -1;
    for (int i=startpos; i < N_; ++i) {
      if (diff_weights_[i] > best_move) {
	best_move = diff_weights_[i];
	best_pos = i;
      }
    }
    if (best_pos < 0 || !ImprovingMove(best_pos)) {
      // No more profitable moves
      break;
    }
    
    // Update the diff_weights_ variable and objective
    UpdateCutValues(best_pos);
  }
}

void ExtendedSolution::AllFirst1Swap(int startpos) {
  // Take all profitable 1-moves, taking the first one we find when scanning
  // the nodes/variables sequentially.
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (int i=startpos; i < N_; ++i) {
      if (ImprovingMove(i)) {
	UpdateCutValues(i);
	move_made = true;
	break;
      }
    }
  }
}

void ExtendedSolution::AllShuffle1Swap(int startpos) {
  // First, build a list of indices to consider and shuffle them
  std::vector<int> indices;
  for (int idx=startpos; idx < N_; ++idx) {
    indices.push_back(idx);
  }
  std::random_shuffle(indices.begin(), indices.end());

  // Take all profitable 1-moves, taking the first one we find when scanning the
  // nodes/variables sequentially.
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (auto iter=indices.begin(); iter != indices.end(); ++iter) {
      if (ImprovingMove(*iter)) {
	UpdateCutValues(*iter);
	move_made = true;
	break;
      }
    }
  }
}

double ExtendedSolution::DiffWeightStandardDeviation() const {
  // Compute the standard deviation in one pass:
  // http://www.strchr.com/standard_deviation_in_one_pass
  double sum = 0.0;
  double sq_sum = 0.0;
  for (int i=0; i < N_; ++i) {
    sum += diff_weights_[i];
    sq_sum += diff_weights_[i] * diff_weights_[i];
  }
  double mean = sum / N_;
  return sqrt(sq_sum / N_ - mean * mean);
}

ExtendedSolution& ExtendedSolution::operator=(const ExtendedSolution &rhs) {
  BaseSolution::operator=(rhs);
  diff_weights_ = rhs.diff_weights_;
  return *this;
}

ExtendedSolution::ExtendedSolution(const ExtendedSolution& x) :
  BaseSolution(x),
  diff_weights_(x.diff_weights_) {}
