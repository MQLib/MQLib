#ifndef HEURISTICS_QUBO_KATAYAMA_2000_H_
#define HEURISTICS_QUBO_KATAYAMA_2000_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Katayama2000QUBOSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Katayama2000QUBOSolution
  Katayama2000QUBOSolution(const QUBOSolution &x) :
   QUBOSolution(x) {}

  static Katayama2000QUBOSolution
    Crossover(const Katayama2000QUBOSolution& x1,
	      const Katayama2000QUBOSolution& x2) {
    return Katayama2000QUBOSolution(x1, x2);
  }

  // Run procedure Variant-k-Opt-Local-Search from Fig. 1
  void VariantKOpt();

  // Flip randomly selected bits in solution
  void Mutate();

 private:
  // Build a new solution by crossover
  Katayama2000QUBOSolution(const Katayama2000QUBOSolution& x1,
			   const Katayama2000QUBOSolution& x2);
};

class Katayama2000Elite {
 public:
  // Build a population of size PS, consisting of random solutions with local
  // search performed on each.
  Katayama2000Elite(const QUBOInstance& qi, int PS, QUBOHeuristic *heuristic);

  // Update the population with a new set of locally optimal child solutions.
  void Update(const std::vector<Katayama2000QUBOSolution>& x);

  // Perform population diversification, if need be.
  void Diversify();

  // getters
  const std::vector<Katayama2000QUBOSolution>& get_P() const  {  return P_; }
  int get_size() const  {  return P_.size(); }
  const Katayama2000QUBOSolution& get_best() const  {  return P_[0]; }

 private:
  // The population of solutions, ordered decreasing by objective value
  std::vector<Katayama2000QUBOSolution> P_;

  int PS_;  // Population size
  QUBOHeuristic *heuristic_;
  int stepsSinceImprovement_;  // How many updates since last improvement?
};

class Katayama2000 : public QUBOHeuristic {
 public:
  Katayama2000(const QUBOInstance& qi, double runtime_limit, bool validation,
               QUBOCallback *qc);
};

#endif
