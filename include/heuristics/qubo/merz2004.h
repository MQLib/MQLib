#ifndef HEURISTICS_QUBO_MERZ_2004_H_
#define HEURISTICS_QUBO_MERZ_2004_H_

#include "heuristics/qubo/merz2002.h"
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

// merz2004 borrows many components from previous papers:
//   - Its randomized greedy procedure is identical to the one from merz2002.
//   - Its mutate procedure is identical to the one from katayama2000, except
//     for a different parameter.
//   - Its randomized k-opt is nearly identical to the one from katayama2000.
//     The only difference is that the k-opt search can be iterated early
//     if there have been a certain number of steps with no improvement.
//   - Its main structure is basically identical to that of katayama2000.
// The only code we'll use directly is the construction procedure, since there
// are changes to the other procedures. Note that this paper has a novel
// crossover procedure, so it's different from previous approaches by these
// authors.

class Merz2004Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Merz2004Solution
 Merz2004Solution(const QUBOSolution& x) :
  QUBOSolution(x) {}
  
  // Generate a new solution via crossover (Fig 7)
  static Merz2004Solution Crossover(const Merz2004Solution& x1,
				    const Merz2004Solution& x2) {
    return Merz2004Solution(x1, x2);
  }

  // Run randomized greedy heuristic from Fig 1 (identical to RandomizedGreedy
  // from Merz2002PartialSolution).
  static Merz2004Solution RandomizedGreedy(const QUBOInstance& qi,
					   QUBOHeuristic *heuristic) {
    return QUBOSolution(Merz2002PartialSolution::RandomizedGreedy(qi,
								  heuristic));
  }

  // Run randomized k-opt from Fig 3 with the limit on the number of inner
  // loop iterations without improvement (identical to VariantKOpt from
  // katayama2000 other than this limitation on the inner loop).
  void RandomizedKOpt();

  // Flip randomly selected bits
  void Mutate();

 private:
  Merz2004Solution(const Merz2004Solution& x1, const Merz2004Solution& x2);
};

// Maintain a set of elite solutions (modeled off Katayama2000Elite)
class Merz2004Elite {
 public:
  // Build a population of size PS, consisting of randomized greedy solutions
  // with local search performed on each.
  Merz2004Elite(const QUBOInstance& qi, int PS, QUBOHeuristic *heuristic);

  // Set P_ to be up to the PS_ best solutions in x, limiting to non-duplicated
  // solutions.
  void SelectNonDuplicated(std::vector<Merz2004Solution>* x);

  // Update the population with a new set of locally optimal child solutions.
  void Update(const std::vector<Merz2004Solution>& x);

  // Perform population diversification, if need be.
  void Diversify();

  // getters
  const std::vector<Merz2004Solution>& get_P() const  {  return P_; }
  const Merz2004Solution& get_best() const  {  return P_[0]; }

 private:
  // The population of solutions, ordered decreasing by objective value
  std::vector<Merz2004Solution> P_;

  int PS_;  // Population size
  int stepsSinceImprovement_;  // How many updates since last improvement?
};

class Merz2004 : public QUBOHeuristic {
 public:
  Merz2004(const QUBOInstance& qi, double runtime_limit, bool validation,
           QUBOCallback *qc);
};

#endif
