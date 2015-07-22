#ifndef HEURISTICS_EXTENDED_SOLUTION_H_
#define HEURISTICS_EXTENDED_SOLUTION_H_

#include <vector>
#include "heuristics/base_solution.h"

class ExtendedSolution : public BaseSolution {
 public:
  // Constructor takes the number of nodes in the graph / vars in the problem
  // and the initial value for the assignments_ vector in the init_assignment
  // parameter.
  ExtendedSolution(int N, int init_assignment);

  // 1-swap functions: startpos means you won't consider 1-swaps for any index
  // less than this parameter.

  // Perform all available 1-moves, at each step selecting the most valuable.
  void AllBest1Swap(int startpos = 0);

  // Perform all 1-moves, at each step selecting the first improving move
  void AllFirst1Swap(int startpos = 0);

  // Perform all 1-moves, iteratively taking the first improving move on a
  // shuffled list of indices.
  void AllShuffle1Swap(int startpos = 0);

  // A popular simulated annealing update mechanism uses the standard deviation
  // of the diff_weights_ values, so we'll provide an accessor for this value.
  double DiffWeightStandardDeviation() const;

  // Identify if a move sufficiently improves the quality of this solution to be
  // taken.
  bool ImprovingMove(int index) const {
    return ImprovesOver(weight_ + diff_weights_[index], weight_);
  }
  bool ImprovingMove(double diff_weight) const {
    return ImprovesOver(weight_ + diff_weight, weight_);
  }
  bool ImprovesOverAfterMove(double weight2, int index) const {
    return ImprovesOver(weight_ + diff_weights_[index], weight2);
  }
  bool ImprovesOverAfterMove(const BaseSolution& other, int index) const {
    return ImprovesOver(weight_ + diff_weights_[index], other.get_weight());
  }
  bool NonDetrimentalMove(int index) const {
    return !ImprovesOver(weight_, weight_ + diff_weights_[index]);
  }

  // Getter
  const std::vector<double>& get_diff_weights() const {  return diff_weights_; }

  // Assignment operator
  ExtendedSolution& operator=(const ExtendedSolution &rhs);

  // Copy constructor
  ExtendedSolution(const ExtendedSolution& x);

 protected:
  // Switch the set of update_index, updating the assignments (x), the
  // diff_weights, and the objective. Extending classes must implement.
  virtual void UpdateCutValues(int update_index, std::vector<int>* x,
			       std::vector<double>* diff_weights,
			       double *objective) const = 0;

  // Switch the set of update_index, updating class variables assignments_,
  // diff_weights_, and weight_
  void UpdateCutValues(int update_index) {
    UpdateCutValues(update_index, &assignments_, &diff_weights_, &weight_);
  }

  // Amount you would gain from switching sets / flipping variable
  std::vector<double> diff_weights_;

 private:
  // Default constructor disabled
  ExtendedSolution();
};

#endif
