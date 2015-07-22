#ifndef HEURISTICS_QUBO_QUBO_SOLUTION_H_
#define HEURISTICS_QUBO_QUBO_SOLUTION_H_

#include <vector>
#include "heuristics/extended_solution.h"
#include "heuristics/qubo/qubo_partial_solution.h"
#include "problem/qubo_heuristic.h"

class QUBOSolution : public ExtendedSolution {
 public:
  // Copy over a QUBOPartialSolution with no fractional assignments.
  QUBOSolution(const QUBOPartialSolution& x);

  // Initialize an all-zero QUBO solution
  static QUBOSolution ZeroSolution(const QUBOInstance& qi,
                                   QUBOHeuristic *heuristic) {
    return QUBOSolution(qi, heuristic, 0, 0);
  }

  // Initialize a completely random solution (p=0.5 for each variable)
  static QUBOSolution RandomSolution(const QUBOInstance& qi,
				     QUBOHeuristic *heuristic) {
    return QUBOSolution(qi, heuristic, 0);
  }

  // Initialize a random solution (p for each variable provided)
  static QUBOSolution RandomSolution(const QUBOInstance& qi,
				     const std::vector<double>& p,
				     QUBOHeuristic *heuristic) {
    return QUBOSolution(qi, p, heuristic);
  }

  // Initialize from vector of assignments
  QUBOSolution(const std::vector<int>& assignments, const QUBOInstance& mi,
               QUBOHeuristic *heuristic);

  // Perform all available 2-moves, at each step selecting the most valuable.
  void AllBest2Swap();

  // Perform all 2-moves, at each step selecting the first improving move
  void AllFirst2Swap();
  
  // Assignment operator
  QUBOSolution& operator=(const QUBOSolution &rhs);

  // Copy constructor
  QUBOSolution(const QUBOSolution& x);

  // This function initializes diff_weights_ and weight_, given that
  // assignments_ is populated.
  void PopulateFromAssignments();

 protected:
  // Initialize qi and heuristic members but nothing else.
  // ********************** IMPORTANT NOTE ***************************
  // diff_weights_ is not properly set after using this constructor, so you
  // must either set that vector manually or by using PopulateFromAssignments
  // *****************************************************************
  QUBOSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic);

  // Switch the set of update_index, updating the assignments (x), the
  // diff_weights, and the objective.
  void UpdateCutValues(int update_index, std::vector<int>* x,
		       std::vector<double>* diff_weights,
		       double *objective) const;

  using ExtendedSolution::UpdateCutValues;  // Unhide single-argument version

  const QUBOInstance& qi_;
  // The associated heuristic (for reporting purposes)
  QUBOHeuristic *heuristic_;
  
 private:
  // Initialize a completely random solution
  QUBOSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic, int ignored);

  // Initialize an all-zero solution
  QUBOSolution(const QUBOInstance& qi, QUBOHeuristic *heuristic, int ignored1,
               int ignored2);

  // Random solution given variable probabilities
  QUBOSolution(const QUBOInstance& qi, const std::vector<double>& p,
	       QUBOHeuristic *heuristic);
};

#endif
