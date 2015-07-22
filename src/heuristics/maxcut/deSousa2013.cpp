#include <iostream>
#include "heuristics/maxcut/deSousa2013.h"
#include "heuristics/maxcut/max_cut_simple_solution.h"
#include "problem/heuristic.h"
#include "util/random.h"

void deSousa2013::AlgorithmIteration() {
  // Parameters (population size and number of generations set based on
  // email communication with paper authors).
  double alpha = 0.5;
  int Ssize = 50;
  int G = 500;

  std::vector<double> p(mi_.get_size(), 0.5);
  MaxCutSimpleSolution best = MaxCutSimpleSolution::RandomSolution(mi_, p, this);
  for (int iter=0; iter < G; ++iter) {
    // Count of assignments in population
    std::vector<int> counts(mi_.get_size(), 0);

    // Update counts from the best solution
    const std::vector<int>& best_assignments = best.get_assignments();
    for (int i=0; i < mi_.get_size(); ++i) {
      if (best_assignments[i] == 1) {
        ++counts[i];
      }
    }

    // Generate the population for this iteration.
    for (int ct=0; ct < Ssize; ++ct) {
      MaxCutSimpleSolution x = MaxCutSimpleSolution::RandomSolution(mi_, p, this);

      // If new solution improves over best, replace the best with it.
      if (x.ImprovesOver(best)) {
        best = x;
        if (!Report(best)) {
          return;
        }
      }

      // Update counts using this solution's assignment
      const std::vector<int>& assignments = x.get_assignments();
      for (int i=0; i < mi_.get_size(); ++i) {
        if (assignments[i] == 1) {
          ++counts[i];
        }
      }
    }

    // Update p based on the counts. If the current p was already basically
    // converged then return (we'll do a random restart).
    bool converged = true;
    for (int i=0; i < mi_.get_size(); ++i) {
      if (p[i] > 0.001 || p[i] < 0.999) {
        converged = false;
      }
      p[i] = (1.0 - alpha) * p[i] + alpha * counts[i] / (Ssize + 1.0);
    }
    if (converged) {
      return;
    }

    // If we are out of runtime, then return.
    if (!Report(best)) {
      return;
    }
  }
}

deSousa2013::deSousa2013(const MaxCutInstance &mi, double runtime_limit,
                         bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  while (MaxCutHeuristic::Report()) {
    AlgorithmIteration();
  }
}
