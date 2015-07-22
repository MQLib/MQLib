#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include "heuristics/base_solution.h"

BaseSolution::BaseSolution(int N, int init_assignment) :
  assignments_(N, init_assignment),
  weight_(0.0),
  N_(N) {}

BaseSolution::BaseSolution(const std::vector<int>& assignments, double weight) :
  assignments_(assignments),
  weight_(weight),
  N_(assignments.size()) {}

int BaseSolution::SymmetricDifference(const BaseSolution& other) const {
  int num_different = 0;
  for (int i=0; i < N_; ++i) {
    num_different += (assignments_[i] != other.assignments_[i]);
  }
  return num_different;
}

int BaseSolution::SymmetricDifference(const BaseSolution& other,
				      std::vector<int>* diff) const {
  diff->clear();
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      diff->push_back(i);
    }
  }
  return diff->size();
}

int BaseSolution::SymmetricDifference(const BaseSolution& other,
				      std::vector<int>* diff,
				      std::vector<int>* common) const {
  diff->clear();
  common->clear();
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      diff->push_back(i);
    } else {
      common->push_back(i);
    }
  }
  return diff->size();
}

bool BaseSolution::operator==(const BaseSolution& other) const {
  // If the weights are different, no need to check the assignments
  if (ImprovesOver(other) || other.ImprovesOver(*this)) {
    return false;
  }

  // Check the assignments
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      return false;
    }
  }
  return true;
}

BaseSolution& BaseSolution::operator=(const BaseSolution &rhs) {
  assignments_ = rhs.assignments_;
  weight_ = rhs.weight_;
  return *this;
}

BaseSolution::BaseSolution(const BaseSolution& x) :
  assignments_(x.assignments_),
  weight_(x.weight_),
  N_(x.N_) {}

bool BaseSolution::ImprovesOver(double weight1, double weight2) {
  return weight1 - weight2 > 1e-6 &&
    (weight2 <= 0.0 || (weight1 / weight2 - 1.0) >= 1e-10);
}

void BaseSolution::PrintSolution() const {
  for (int i=0; i < N_; ++i) {
    if (i != 0) {
      std::cout << " ";
    }
    std::cout << assignments_[i];
  }
  std::cout << std::endl;
}
