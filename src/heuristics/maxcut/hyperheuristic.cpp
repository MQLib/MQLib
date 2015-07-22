#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "heuristics/heuristic_factory.h"
#include "heuristics/maxcut/hyperheuristic.h"
#include "metrics/max_cut_metrics.h"
#include "util/random.h"
#include "util/randomForest.h"

bool MaxCutHyperheuristic::FileExists(const std::string& filename) {
  FILE *f = fopen(filename.c_str(), "r");
  if (f) {
    fclose(f);
    return true;
  } else {
    return false;
  }
}

void MaxCutHyperheuristic::UpdateBestModel(std::string code, Prob problem,
                                           const std::vector<double>& metrics,
                                           double* bestProbability,
                                           Prob* bestProblem,
                                           std::string* bestCode,
                                           int* numBest) {
  std::ostringstream fname;
  fname << "hhdata/" << code << ".rf";
  std::string filename = fname.str();
  if (FileExists(filename)) {
    RandomForest rf(filename);
    double probability = rf.Predict(metrics);
    if (probability > *bestProbability) {
      // New best
      *bestProbability = probability;
      *bestProblem = problem;
      *bestCode = code;
      *numBest = 1;
    } else if (probability == *bestProbability &&
               Random::RandInt(0, *numBest) == *numBest) {
      // Tied the best and selected by streaming algorithm
      *bestProblem = problem;
      *bestCode = code;
      ++(*numBest);
    }
  }
}

HyperheuristicMaxCutCallback::HyperheuristicMaxCutCallback(MaxCutHyperheuristic* mch) :
  mch_(mch) {}

bool HyperheuristicMaxCutCallback::Report(const MaxCutSimpleSolution& solution,
                                          bool newBest, double runtime) {
  return mch_->Report(solution);
}

bool HyperheuristicMaxCutCallback::Report(const MaxCutSimpleSolution& solution,
                                          bool newBest, double runtime,
                                          int iter) {
  return mch_->Report(solution, iter);
}

HyperheuristicQUBOCallback::HyperheuristicQUBOCallback(MaxCutHyperheuristic* mch, const MaxCutInstance& mi) :
  mch_(mch),
  mi_(mi) {}

bool HyperheuristicQUBOCallback::Report(const QUBOSimpleSolution& solution,
                                        bool newBest, double runtime) {
  if (newBest) {
    // Need to convert from QUBO to Max-Cut and report
    MaxCutSimpleSolution mcSol(solution, mi_, mch_);
    return mch_->Report(mcSol);
  } else {
    return mch_->Report();  // Not new best solution, so just check term. crit.
  }
}

bool HyperheuristicQUBOCallback::Report(const QUBOSimpleSolution& solution,
                                        bool newBest, double runtime,
                                        int iter) {
  if (newBest) {
    // Need to convert from QUBO to Max-Cut and report
    MaxCutSimpleSolution mcSol(solution, mi_, mch_);
    return mch_->Report(mcSol, iter);
  } else {
    return mch_->Report(iter);  // Not new best solution, so check term. crit.
  }
}

MaxCutHyperheuristic::MaxCutHyperheuristic(const MaxCutInstance&mi,
                                           double runtime_limit,
                                           bool validation,
                                           MaxCutCallback *mc, int seed,
                                           std::string* selected) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Step 1: Calculate graph metrics for this instance.
  GraphMetrics gm(mi);
  std::vector<double> metrics;
  gm.AllMetrics(&metrics, NULL);

  // Step 2: Obtain predicted probabilities from each random forest model.
  // Best-performing model information (using streaming alg to select best)
  double bestProbability = -1.0;
  Prob bestProblem = MaxCut;
  std::string bestCode = "";
  int numBest = 1;  // Number tied for best

  // Check the Max-Cut heuristics
  HeuristicFactory factory;
  std::vector<std::string> codes;
  factory.MaxCutHeuristicCodes(&codes);
  for (int i=0; i < codes.size(); ++i) {
    UpdateBestModel(codes[i], MaxCut, metrics, &bestProbability, &bestProblem,
                    &bestCode, &numBest);
  }

  // Check the QUBO heuristics
  factory.QUBOHeuristicCodes(&codes);
  for (int i=0; i < codes.size(); ++i) {
    UpdateBestModel(codes[i], QUBO, metrics, &bestProbability, &bestProblem,
                    &bestCode, &numBest);
  }
  if (selected) {
    *selected = bestCode;
  }

  // Step 3: Run the selected heuristic "H", using callbacks to capture all the
  // reported solutions from "H" and report them for the hyper-heuristic.
  // Because several previous steps may have used random draws or set the
  // random seed, re-set the seed to the original here.
  srand(seed);
  if (bestProblem == MaxCut) {
    // Using a Max-Cut heuristic
    HyperheuristicMaxCutCallback callback(this);
    // Run with our callback and no validation (solutions will be validated
    // with the hyperheuristic, so no need to double validate)
    Heuristic *h = factory.RunMaxCutHeuristic(bestCode, mi, runtime_limit,
                                              false, &callback);
    delete h;  // We don't need to keep around the pointer
  } else if (bestProblem == QUBO) {
    // Using a QUBO heuristic
    HyperheuristicQUBOCallback callback(this, mi);
    QUBOInstance qi(mi);  // Need to construct a QUBOInstance for heuristic
    // Run with our callback and no validation (solutions will be validated
    // with the hyperheuristic, so no need to double validate)
    Heuristic *h = factory.RunQUBOHeuristic(bestCode, qi, runtime_limit,
                                            false, &callback);
    delete h;  // We don't need to keep around the pointer
  }
}
