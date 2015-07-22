#ifndef HEURISTICS_MAXCUT_DUARTE_2005_H_
#define HEURISTICS_MAXCUT_DUARTE_2005_H_

#include "problem/max_cut_heuristic.h"
#include "heuristics/maxcut/max_cut_solution.h"

class Duarte2005Solution : public MaxCutSolution {
 public:
  // Convert from MaxCutSolution to Duarte2005Solution
  Duarte2005Solution(const MaxCutSolution &x);

  // Genetic crossover operator
  static Duarte2005Solution FixCross(const Duarte2005Solution& father,
                                     const Duarte2005Solution& mother) {
    return Duarte2005Solution(father, mother);
  }

  // Greedy local search
  void Greedy1Swap();

  // Mutation
  void Mutate(double flip_chance);

  // VNS that uses Greedy1Swap
  void VNS(int k_max);

 private:
  // Genetic crossover operator constructor
  Duarte2005Solution(const Duarte2005Solution& father,
                     const Duarte2005Solution& mother);
};


class Duarte2005 : public MaxCutHeuristic {
 public:
  Duarte2005(const MaxCutInstance& mi, double runtime_limit, bool validation,
             MaxCutCallback *mc);
};

#endif
