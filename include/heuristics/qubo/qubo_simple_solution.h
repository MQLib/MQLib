#ifndef HEURISTICS_QUBO_QUBO_SIMPLE_SOLUTION_H_
#define HEURISTICS_QUBO_QUBO_SIMPLE_SOLUTION_H_

#include <vector>
#include "heuristics/base_solution.h"
class QUBOHeuristic;
class QUBOInstance;
#include "problem/qubo_heuristic.h"
class MaxCutSimpleSolution;
#include "heuristics/maxcut/max_cut_simple_solution.h"

class QUBOSimpleSolution : public BaseSolution {
 public:
  // Initialize a completely random solution (p=0.5 for each vertex)
  static QUBOSimpleSolution RandomSolution(const QUBOInstance& qi,
                                           QUBOHeuristic *heuristic) {
    return QUBOSimpleSolution(qi, heuristic, 0, 0);  // Private constructor
  }

  // Initialize a random solution (p for each vertex provided)
  static QUBOSimpleSolution RandomSolution(const QUBOInstance& qi,
                                           const std::vector<double>& p,
                                           QUBOHeuristic *heuristic) {
    return QUBOSimpleSolution(qi, p, heuristic);
  }

  // Initialize qi and heuristic members but nothing else
  QUBOSimpleSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic,
                     int initVal = 0);

  // Initialize from qi, heuristic, assignments, and weight
  QUBOSimpleSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic,
                     const std::vector<int>& assignments, double weight);

  // Initialize from a Max-Cut solution sol and QUBOInstance qi, under the
  // assumption that sol was obtained by solving a MaxCutInstance mi that was
  // constructed by reducing qi to mi. This is a light-weight constructor that
  // does not re-compute the objective value, merely applying the problem
  // reduction to convert the instance.
  QUBOSimpleSolution(const MaxCutSimpleSolution& sol, const QUBOInstance& qi,
                     QUBOHeuristic *heuristic);

  // Assignment operator
  QUBOSimpleSolution& operator=(const QUBOSimpleSolution &rhs);

  // Copy constructor
  QUBOSimpleSolution(const QUBOSimpleSolution& x);

  // This function initializes weight_, given that assignments_ is populated.
  void PopulateFromAssignments();

 protected:
  const QUBOInstance& qi_;
  // The associated heuristic (for reporting purposes)
  QUBOHeuristic *heuristic_;

 private:
  // Initialize a completely random solution
  QUBOSimpleSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic,
                     int ignored1, int ignored2);

  // Random solution given vertex probabilities
  QUBOSimpleSolution(const QUBOInstance& qi, const std::vector<double>& p,
                     QUBOHeuristic *heuristic);
};

#endif
