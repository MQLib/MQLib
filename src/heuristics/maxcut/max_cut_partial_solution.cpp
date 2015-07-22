#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "heuristics/maxcut/max_cut_partial_solution.h"
#include "problem/max_cut_instance.h"

MaxCutPartialSolution::MaxCutPartialSolution(const MaxCutInstance& mi,
					     MaxCutHeuristic *heuristic) :
  mi_(mi),
  heuristic_(heuristic),
  N_(mi.get_size()),
  assignments_(N_, 0),
  gainS_(N_, 0.0),
  gainNS_(N_, 0.0),
  num_unassigned_(N_),
  weight_(0.0) {}

void MaxCutPartialSolution::PopulateFromAssignments() {
  // Initialize gainS_, gainNS_, num_unassigned_ and weight_
  gainS_.assign(N_, 0.0);
  gainNS_.assign(N_, 0.0);
  num_unassigned_ = 0;
  weight_ = 0.0;

  // Set num_unassigned_ to proper value
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] == 0) {
      ++num_unassigned_;
    }
  }

  // Iterate through the edges, computing gainS_, gainNS_ and weight_.
  for (auto iter = mi_.get_all_edges_begin();
       iter != mi_.get_all_edges_end(); ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double wij = iter->second;
    switch (assignments_[i]) {
      // i is in set NS:
    case -1:
      switch (assignments_[j]) {
      case -1:
	gainS_[i] += wij;
	gainS_[j] += wij;
	break;
      case 0:
	gainS_[j] += wij;
	break;
      case 1:
	weight_ += wij;
	gainS_[i] -= wij;
	gainNS_[j] -= wij;
	break;
      default:
	std::cout << "Illegal assignment:" << assignments_[j] << std::endl;
	exit(0);
      }
      break;

      // i is not assigned to a set:
    case 0:
      switch (assignments_[j]) {
      case -1:
	gainS_[i] += wij;
	break;
      case 0:
	break;
      case 1:
	gainNS_[i] += wij;
	break;
      default:
	std::cout << "Illegal assignment:" << assignments_[j] << std::endl;
	exit(0);
      }
      break;

      // i is in set S:
    case 1:
      switch (assignments_[j]) {
      case -1:
	weight_ += wij;
	gainNS_[i] -= wij;
	gainS_[j] -= wij;
	break;
      case 0:
	gainNS_[j] += wij;
	break;
      case 1:
	gainNS_[i] += wij;
	gainNS_[j] += wij;
	break;
      default:
	std::cout << "Illegal assignment:" << assignments_[j] << std::endl;
	exit(0);
      }
      break;
    default:
      std::cout << "Illegal assignment:" << assignments_[i] << std::endl;
      exit(0);
    }
  }
}

void MaxCutPartialSolution::UpdateCutValues(int update_index, int new_value) {
  if ((new_value != -1 && new_value != 1) || update_index < 0 ||
      update_index >= N_) {
    std::cout << "Illegal parameters to UpdateCutValues" << std::endl;
    exit(0);
  }
  
  // If this isn't actually an update, just return
  if (new_value == assignments_[update_index]) {
    return;
  }

  // If the value at update_index was 0 (unassigned), then we're decreasing the
  // number of unassigned values by 1 here.
  if (assignments_[update_index] == 0) {
    --num_unassigned_;
  }

  // Iterate the edges incident to update_index, updating gainS_ and gainNS_ for
  // the incident nodes. We'll loop separately for each type of update to
  // update_index (0 -> -1, 0 -> 1, -1 -> 1, and 1 -> -1).
  if (assignments_[update_index] == 0 && new_value == -1) {
    for (auto iter = mi_.get_edges_begin(update_index);
	 iter != mi_.get_edges_end(update_index); ++iter) {
      int j = iter->first;
      double wij = iter->second;
      if (assignments_[j] == 1) {
	gainNS_[j] -= wij;
      } else {
	gainS_[j] += wij;
      }
    }
  } else if (assignments_[update_index] == 0 && new_value == 1) {
    for (auto iter = mi_.get_edges_begin(update_index);
	 iter != mi_.get_edges_end(update_index); ++iter) {
      int j = iter->first;
      double wij = iter->second;
      if (assignments_[j] == -1) {
	gainS_[j] -= wij;
      } else {
	gainNS_[j] += wij;
      }
    }
  } else if (assignments_[update_index] == -1 && new_value == 1) {
    for (auto iter = mi_.get_edges_begin(update_index);
	 iter != mi_.get_edges_end(update_index); ++iter) {
      int j = iter->first;
      double wij = iter->second;

      switch (assignments_[j]) {
      case -1:
	gainS_[j] -= 2.0 * wij;
	break;
      case 0:
	gainS_[j] -= wij;
	gainNS_[j] += wij;
	break;
      default:  // j == 1
	gainNS_[j] += 2.0 * wij;
      }
    }
  } else {
    for (auto iter = mi_.get_edges_begin(update_index);
	 iter != mi_.get_edges_end(update_index); ++iter) {
      int j = iter->first;
      double wij = iter->second;

      switch (assignments_[j]) {
      case -1:
	gainS_[j] += 2.0 * wij;
	break;
      case 0:
	gainS_[j] += wij;
	gainNS_[j] -= wij;
	break;
      default:  // j == 1
	gainNS_[j] -= 2.0 * wij;
      }
    }
  }

  // Update weight_ as well as assignments_, gainS_ and gainNS_ for update_index
  if (new_value == -1) {
    assignments_[update_index] = -1;
    weight_ += gainNS_[update_index];
    gainS_[update_index] -= gainNS_[update_index];
    gainNS_[update_index] = 0.0;
  } else {
    assignments_[update_index] = 1;
    weight_ += gainS_[update_index];
    gainNS_[update_index] -= gainS_[update_index];
    gainS_[update_index] = 0.0;
  }
}
