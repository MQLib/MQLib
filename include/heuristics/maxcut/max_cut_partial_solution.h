#ifndef HEURISTICS_MAXCUT_MAXCUT_PARTIAL_SOLUTION_H_
#define HEURISTICS_MAXCUT_MAXCUT_PARTIAL_SOLUTION_H_

#include <vector>
#include "problem/max_cut_heuristic.h"

// Some constructive procedures for MAXCUT assign edges one at a time
// and maintain the gains associated with assigning a node either to S
// (assignment value 1) or NS (assignment value -1). The MaxCutPartialSolution
// supports such procedures, handling variable updates. The MaxCutPartialSolution
// can be converted (cheaply) to a MaxCutSolution through a constructor in
// MaxCutSolution.
class MaxCutPartialSolution {
 public:
  // Getters
  double get_weight() const {  return weight_; }
  int get_num_unassigned() const {  return num_unassigned_; }
  const std::vector<int>& get_assignments() const {  return assignments_; }
  const std::vector<double>& get_gainS() const {  return gainS_; }
  const std::vector<double>& get_gainNS() const {  return gainNS_; }
  const MaxCutInstance& get_mi() const {  return mi_; }
  MaxCutHeuristic* get_heuristic() const {  return heuristic_; }

 protected:
  // Initialize a solution to all unassigned
  MaxCutPartialSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic);

  // Update gainS_, gainNS_, num_unassigned_, and weight_ from assignments_
  void PopulateFromAssignments();

  // Set the value of update_index to the specified new_value, updating
  // gainS_, gainNS_, assignments_, num_unassigned_, and weight_.
  void UpdateCutValues(int update_index, int new_value);

  // Problem instance
  const MaxCutInstance& mi_;
  // Associated heuristic
  MaxCutHeuristic *heuristic_;
  // Convenience variable for number of variables in the problem
  int N_;
  // Assignment of each node (1 means S, -1 means NS, 0 means unassigned).
  std::vector<int> assignments_;
  // Change in objective from assigning each vertex to set S (value 1)
  std::vector<double> gainS_;
  // Change in objective from assigning each vertex to set NS (value -1)
  std::vector<double> gainNS_;
  // Number of assignments that are unassigned (value 0)
  int num_unassigned_;
  // Objective value
  double weight_;
};

#endif
