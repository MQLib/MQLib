#ifndef HEURISTICS_MAXCUT_MAX_CUT_SIMPLE_SOLUTION_H_
#define HEURISTICS_MAXCUT_MAX_CUT_SIMPLE_SOLUTION_H_

#include <vector>
#include "heuristics/base_solution.h"
class MaxCutInstance;
class MaxCutHeuristic;
#include "problem/max_cut_heuristic.h"
class QUBOSimpleSolution;
#include "heuristics/qubo/qubo_simple_solution.h"

class MaxCutSimpleSolution : public BaseSolution {
 public:
  // Initialize a completely random solution (p=0.5 for each vertex)
  static MaxCutSimpleSolution RandomSolution(const MaxCutInstance& mi,
                                             MaxCutHeuristic *heuristic) {
    return MaxCutSimpleSolution(mi, heuristic, 0, 0);  // Private constructor
  }

  // Initialize a random solution (p for each vertex provided)
  static MaxCutSimpleSolution RandomSolution(const MaxCutInstance& mi,
                                             const std::vector<double>& p,
                                             MaxCutHeuristic *heuristic) {
    return MaxCutSimpleSolution(mi, p, heuristic);
  }

  // Initialize mi and heuristic members but nothing else. This is a light-
  // weight constructor that labels the solution as having weight 0.
  MaxCutSimpleSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
                       int initVal = -1);

  // Initialize from mi, heuristic, assignments, and weight. This is a light-
  // weight constructor that does not re-compute the weight.
  MaxCutSimpleSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
                       const std::vector<int>& assignments, double weight);

  // Initialize from a QUBO solution sol and MaxCutInstance mi, under the
  // assumption that sol was obtained by solving a QUBOInstance qi that was
  // constructed by reducing mi to qi. This is a light-weight constructor that
  // does not re-compute the objective value, merely applying the problem
  // reduction to convert the instance.
  MaxCutSimpleSolution(const QUBOSimpleSolution& sol, const MaxCutInstance& mi,
                       MaxCutHeuristic *heuristic);

  // Assignment operator
  MaxCutSimpleSolution& operator=(const MaxCutSimpleSolution &rhs);

  // Copy constructor
  MaxCutSimpleSolution(const MaxCutSimpleSolution& x);

  // This function initializes weight_, given that assignments_ is populated.
  void PopulateFromAssignments();

 protected:
  const MaxCutInstance& mi_;
  // The associated heuristic (for reporting purposes)
  MaxCutHeuristic *heuristic_;

 private:
  // Initialize a completely random solution
  MaxCutSimpleSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
                       int ignored1, int ignored2);

  // Random solution given vertex probabilities
  MaxCutSimpleSolution(const MaxCutInstance& mi, const std::vector<double>& p,
                       MaxCutHeuristic *heuristic);
};

#endif
