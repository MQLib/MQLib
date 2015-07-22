#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "heuristics/qubo/qubo_partial_solution.h"
#include "problem/qubo_instance.h"

QUBOPartialSolution::QUBOPartialSolution(const QUBOInstance& qi,
					 QUBOHeuristic *heuristic) :
  qi_(qi),
  heuristic_(heuristic),
  N_(qi.get_size()),
  assignments_(N_, 0.0),
  diff0_(N_, 0.0),
  diff1_(qi.get_lin()),
  num_frac_(0),
  weight_(0.0) {}

void QUBOPartialSolution::PopulateFromAssignments() {
  // Initialize diff0_, diff1_, num_frac_ and weight_
  diff0_.assign(N_, 0.0);
  diff1_.assign(N_, 0.0);
  num_frac_ = 0;
  weight_ = 0.0;

  // Handle linear terms for updates to weight_, diff0_ and diff1_. Also update
  // num_frac_.
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] < 0.0 || assignments_[i] > 1.0) {
      std::cout << "Invalid assignment in QUBOPartialSolution" << std::endl;
      exit(0);
    } else if (assignments_[i] == 0.0) {
      diff1_[i] += qi_.get_lin()[i];
    } else if (assignments_[i] == 1.0) {
      diff0_[i] -= qi_.get_lin()[i];
      weight_ += qi_.get_lin()[i];
    } else {
      ++num_frac_;
      double prod = assignments_[i] * assignments_[i];
      diff0_[i] -= prod * qi_.get_lin()[i];
      diff1_[i] += (1.0 - prod) * qi_.get_lin()[i];
      weight_ += prod * qi_.get_lin()[i];
    }
  }

  // Iterate through the edges, completing diff0_, diff1_ and weight_.
  for (auto iter = qi_.get_all_nonzero_begin();
       iter != qi_.get_all_nonzero_end(); ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double wij = iter->second;
    double prod = assignments_[i] * assignments_[j];
    diff0_[i] -= prod * wij * 2.0;
    diff0_[j] -= prod * wij * 2.0;
    diff1_[i] += (assignments_[j] - prod) * wij * 2.0;
    diff1_[j] += (assignments_[i] - prod) * wij * 2.0;
    weight_ += prod * wij * 2.0;
  }  
}

void QUBOPartialSolution::UpdateCutValues(int update_index, int new_value) {
  if ((new_value != 0 && new_value != 1) || update_index < 0 ||
      update_index >= N_) {
    std::cout << "Illegal parameters to UpdateCutValues" << std::endl;
    exit(0);
  }
  
  // If the value at update_index was fractional, then we're decreasing the
  // number of fractional indices by 1 here.
  if (assignments_[update_index] != 0.0 && assignments_[update_index] != 1.0) {
    --num_frac_;
  }

  // Iterate the edges incident to update_index, updating diff0_ and diff1_ for
  // the incident nodes
  for (auto iter = qi_.get_nonzero_begin(update_index);
       iter != qi_.get_nonzero_end(update_index); ++iter) {
    int j = iter->first;
    double qij = iter->second;
    if (new_value == 0) {
      diff0_[j] += 2.0 * assignments_[update_index] * assignments_[j] * qij;
      diff1_[j] -=
	2.0 * assignments_[update_index] * (1.0 - assignments_[j]) * qij;
    } else {
      diff0_[j] +=
	2.0 * assignments_[j] * qij * (assignments_[update_index] - 1.0);
      diff1_[j] +=
	2.0*qij * (1.0 - assignments_[j]) * (1.0 - assignments_[update_index]);
    }
  }

  // Update weight_ as well as assignments_, diff0_ and diff1_ for update_index
  if (new_value == 0) {
    assignments_[update_index] = 0.0;
    weight_ += diff0_[update_index];
    diff1_[update_index] -= diff0_[update_index];
    diff0_[update_index] = 0.0;
  } else {
    assignments_[update_index] = 1.0;
    weight_ += diff1_[update_index];
    diff0_[update_index] -= diff1_[update_index];
    diff1_[update_index] = 0.0;
  }
}
