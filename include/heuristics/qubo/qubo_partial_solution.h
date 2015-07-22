#ifndef HEURISTICS_QUBO_QUBO_PARTIAL_SOLUTION_H_
#define HEURISTICS_QUBO_QUBO_PARTIAL_SOLUTION_H_

#include <vector>
#include "problem/qubo_heuristic.h"

// Some constructive procedures for QUBO assign fractional values to variables
// and maintain the gains associated with a switch to either value 0 or value 1.
// QUBOPartialSolution supports such procedures, handling variable updates. The
// QUBOPartialSolution can be converted (cheaply) to a QUBOSolution through a
// constructor in QUBOSolution.
class QUBOPartialSolution {
 public:
  // Getters
  double get_weight() const {  return weight_; }
  int get_num_frac() const {  return num_frac_; }
  const std::vector<double>& get_assignments() const {  return assignments_; }
  const std::vector<double>& get_diff0() const {  return diff0_; }
  const std::vector<double>& get_diff1() const {  return diff1_; }
  const QUBOInstance& get_qi() const {  return qi_; }
  QUBOHeuristic* get_heuristic() const {  return heuristic_; }

 protected:
  // Initialize a solution to all 0 assignments
  QUBOPartialSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic);

  // Update diff0_, diff1_, num_frac_, and weight_ from assignments_
  void PopulateFromAssignments();

  // Set the value of update_index to the specified new_value, updating
  // diff0_, diff1_, assignments_, num_frac_, and weight_.
  void UpdateCutValues(int update_index, int new_value);

  // Problem instance
  const QUBOInstance& qi_;
  // Associated heuristic
  QUBOHeuristic *heuristic_;
  // Convenience variable for number of variables in the problem
  int N_;
  // Assignment of each node (any value between 0 and 1 is valid).
  std::vector<double> assignments_;
  // Change in objective from changing each index to value 0
  std::vector<double> diff0_;
  // Change in objective from changing each index to value 1
  std::vector<double> diff1_;
  // Number of assignments that are fractional (not 0 or 1)
  int num_frac_;
  // Objective value
  double weight_;
};

#endif
