#ifndef HEURISTICS_MAXCUT_HYPERHEURISTIC_H_
#define HEURISTICS_MAXCUT_HYPERHEURISTIC_H_

#include <string>
#include "problem/max_cut_heuristic.h"
#include "problem/qubo_heuristic.h"

class MaxCutHyperheuristic : public MaxCutHeuristic {
 public:
  MaxCutHyperheuristic(const MaxCutInstance&mi, double runtime_limit,
                       bool validation, MaxCutCallback *mc, int seed,
                       std::string* selected);

 private:
  enum Prob {
    MaxCut,
    QUBO
  };

  // Does the indicated file name point to an accessible file?
  bool FileExists(const std::string& filename);

  // Determine if the associated random forest model outperforms the best
  // identified so far, using a streaming algorithm to break ties randomly.
  void UpdateBestModel(std::string code, Prob problem,
                       const std::vector<double>& metrics,
                       double* bestProbability, Prob* bestProblem,
                       std::string* bestCode, int* numBest);
};

// When the hyper-heuristic selects a Max-Cut heuristic, we use this callback to
// capture all reported solutions and appropriately report them for the
// hyper-heuristic.
class HyperheuristicMaxCutCallback : public MaxCutCallback {
 public:
  HyperheuristicMaxCutCallback(MaxCutHyperheuristic* mch);
  bool Report(const MaxCutSimpleSolution& solution, bool newBest,
              double runtime);
  bool Report(const MaxCutSimpleSolution& solution, bool newBest,
              double runtime, int iter);

 private:
  // Pointer to hyperheuristic (used to report solutions on behalf of the
  // hyperheuristic)
  MaxCutHyperheuristic* mch_;
};

// When the hyper-heuristic selects a QUBO heuristic, we use this callback to
// capture all reported solutions and appropriately report them for the
// hyper-heuristic. Whenever a new best solution is found we convert it to a
// Max-Cut solution and report it; otherwise we report without a new solution.
class HyperheuristicQUBOCallback : public QUBOCallback {
 public:
  HyperheuristicQUBOCallback(MaxCutHyperheuristic* mch,
                             const MaxCutInstance& mi);
  bool Report(const QUBOSimpleSolution& solution, bool newBest, double runtime);
  bool Report(const QUBOSimpleSolution& solution, bool newBest,
              double runtime, int iter);

 private:
  // Pointer to hyperheuristic (used to report solutions on behalf of the
  // hyperheuristic)
  MaxCutHyperheuristic* mch_;

  // Pointer to the Max-Cut instance for the hyperheuristic (used to convert
  // passed QUBOSimpleSolution objects into MaxCutSimpleSolution objects).
  const MaxCutInstance& mi_;
};

#endif
