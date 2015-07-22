#ifndef HEURISTICS_BASE_SOLUTION_H_
#define HEURISTICS_BASE_SOLUTION_H_

#include <vector>

class BaseSolution {
 public:
  // Constructor takes the number of nodes in the graph / vars in the problem
  // and the initial value for the assignments_ vector in the init_assignment
  // parameter.
  BaseSolution(int N, int init_assignment);
  
  // Constructor takes the assignments, weight, and problem.
  BaseSolution(const std::vector<int>& assignments, double weight);

  // Return the Hamming distance to another solution.
  int SymmetricDifference(const BaseSolution& other) const;
  
  // Return the Hamming distance to another solution, and populate a vector of
  // the indices that differ.
  int SymmetricDifference(const BaseSolution& other,
			  std::vector<int>* diff) const;

  // Return the Hamming distance to another solution, and populate a vector of
  // the indices that differ and indices in common
  int SymmetricDifference(const BaseSolution& other,
			  std::vector<int>* diff,
			  std::vector<int>* common) const;

  // Identify if one solution improves over another solution. These functions take
  // into account the fact that our problems are stated with floating point
  // numbers instead of integers as weights.
  static bool ImprovesOver(double weight1, double weight2);
  bool ImprovesOver(double other_weight) const {
    return ImprovesOver(weight_, other_weight);
  }
  bool ImprovesOver(const BaseSolution& other) const {
    return ImprovesOver(weight_, other.weight_);
  }

  // Print out the decision variables for a solutio
  void PrintSolution() const;

  // Getters
  double get_weight() const {  return weight_; }
  const std::vector<int>& get_assignments() const {  return assignments_; }

  // Assignment operator
  BaseSolution& operator=(const BaseSolution &rhs);

  // Copy constructor
  BaseSolution(const BaseSolution& x);

  // Equals operator (checks if the solutions have identical assignments_)
  bool operator==(const BaseSolution& other) const;

  // Not equals operator
  bool operator!=(const BaseSolution& other) const {  return !(*this == other); }

  // >, <, >=, <=: compare based on weight_
  bool operator>(const BaseSolution& other) const {
    return weight_ > other.weight_;
  }
  bool operator<(const BaseSolution& other) const {
    return weight_ < other.weight_;
  }
  bool operator>=(const BaseSolution& other) const {
    return weight_ >= other.weight_;
  }
  bool operator<=(const BaseSolution& other) const {
    return weight_ <= other.weight_;
  }

 protected:
  // Assignment of each node (-1/1 in MAXCUT, 0/1 in QUBO)
  std::vector<int> assignments_;
  // Objective value
  double weight_;
  // The number of nodes in the graph / variables in the problem
  int N_;

 private:
  // Default constructor disabled
  BaseSolution();
};

#endif
