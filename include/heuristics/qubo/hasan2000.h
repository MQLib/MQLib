#ifndef HEURISTICS_QUBO_HASAN_2000_H_
#define HEURISTICS_QUBO_HASAN_2000_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

// hasan2000 has a simulated annealing method that is identical to the use
// presented in alkhamis1998, except that a few parameter values are different.
// As a result, we will only implement the genetic algorithm and tabu search.

class Hasan2000Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Hasan2000Solution
 Hasan2000Solution(const QUBOSolution& x) :
  QUBOSolution(x) {}

  // Generate a new solution via one-point crossover with a 2-bit crossing
  // length (Sec. 4.1.3)
  static Hasan2000Solution Crossover(const Hasan2000Solution& x1,
				     const Hasan2000Solution& x2) {
    return Hasan2000Solution(x1, x2);
  }

  // Perform mutation a solution (Sec. 4.1.4)
  void Mutate();

  // Perform the tabu search described in Sec. 5
  void TS();

 private:
  Hasan2000Solution(const Hasan2000Solution& x1, const Hasan2000Solution& x2);
};

// Maintain a set of elite solutions
class Hasan2000Elite {
 public:
  // Build a population of size POP, consisting of random solutions with local
  // search performed on each. Report each of these initial solutions.
  Hasan2000Elite(const QUBOInstance& qi, int POP, QUBOHeuristic *heuristic);

  // Update the population with a new locally optimal child solution, replacing
  // the worst. Only add the solution if it's not already a member of the elite
  // set. Returns if we added the solution.
  bool Update(const Hasan2000Solution& x);

  // Getters
  const std::vector<Hasan2000Solution>& get_P() const  {  return P_; }

 private:
  // The population of solutions.
  std::vector<Hasan2000Solution> P_;

  int POP_;  // Population size
};

class Hasan2000GA : public QUBOHeuristic {
 public:
  Hasan2000GA(const QUBOInstance& qi, double runtime_limit, bool validation,
              QUBOCallback *qc);
};

class Hasan2000TS : public QUBOHeuristic {
 public:
  Hasan2000TS(const QUBOInstance& qi, double runtime_limit, bool validation,
              QUBOCallback *qc);
};

#endif
