#ifndef HEURISTICS_QUBO_LU_2010_H_
#define HEURISTICS_QUBO_LU_2010_H_

#include "heuristics/qubo/qubo_partial_solution.h"
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"
#include <vector>

class Lu2010QUBOSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Lu2010QUBOSolution
  Lu2010QUBOSolution(const QUBOSolution &x);

  // Runs Combination_Operator from Sec 2.5 (either uniform or DG/PR)
  static Lu2010QUBOSolution Combine(const QUBOInstance& qi,
				    const Lu2010QUBOSolution& xj,
				    const Lu2010QUBOSolution& xk,
				    QUBOHeuristic *heuristic) {
    return Lu2010QUBOSolution(qi, xj, xk, heuristic);
  }

  // Perform a tabu search, updating solution to the best TS value
  void TabuSearch();

 private:
  Lu2010QUBOSolution(const QUBOInstance& qi, const Lu2010QUBOSolution& xj,
		     const Lu2010QUBOSolution& xk, QUBOHeuristic *heuristic);
};

// Lu2010PartialSolution enables the DG/PR, which involves variables set to
// fractional values.
class Lu2010PartialSolution : public QUBOPartialSolution {
 public:
  // Generate a solution via DG/PR. By the end of this constructor, there are
  // no fractional variables.
  static Lu2010PartialSolution DGPR(const QUBOInstance& qi,
				    const Lu2010QUBOSolution& xj,
				    const Lu2010QUBOSolution& xk,
				    QUBOHeuristic *heuristic) {
    return Lu2010PartialSolution(qi, xj, xk, heuristic);
  }

 private:
  Lu2010PartialSolution(const QUBOInstance& qi, const Lu2010QUBOSolution& xj,
			const Lu2010QUBOSolution& xk, QUBOHeuristic *heuristic);
};

// Maintain a population of solutions, including the hamming distances between
// each pair of solutions in the population.
class Lu2010Population {
 public:
  // Load up an initial population of size p (each element is random, with tabu
  // search applied)
  Lu2010Population(int p, const QUBOInstance& qi, QUBOHeuristic *heuristic);

  // Section 2.5.1: Randomly select parents with above average hamming distance
  std::pair<const Lu2010QUBOSolution&, const Lu2010QUBOSolution&>
    RandomParents() const;

  bool isConstant() const  {  return avg_HD_ == 0.0;  }

  // Section 2.6: Pool_Updating method
  void UpdatePool(const Lu2010QUBOSolution& x0);

  Lu2010QUBOSolution* GetSolution(int idx)  {  return &(P_[idx]);  }

 private:
  // The population of all solutions
  std::vector<Lu2010QUBOSolution> P_;

  // The population size
  int p_;

  // A p x p matrix of Hamming distances; dist(i, j) is stored in position
  // i*p + j.
  std::vector<int> HD_;

  // The average hamming distance between pairs of solutions
  double avg_HD_;

  // A p x p matrix of distances between solutions in the population (based
  // on variable importance). dist(i, j) is stored in position i*p + j
  std::vector<double> NHD_;

  // The minimum non-Hamming distance (the variable importance one) from each
  // member of the population to all the other ones.
  std::vector<double> min_NHD_;

  // Variable importances
  std::vector<double> VI_;
};


class Lu2010 : public QUBOHeuristic {
 public:
  Lu2010(const QUBOInstance& qi, double runtime_limit, bool validation,
         QUBOCallback *qc);
};

#endif
