#include <math.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "heuristics/maxcut/laguna2009.h"
#include "util/random.h"

Laguna2009CE::Laguna2009CE(const MaxCutInstance& mi, double runtime_limit,
			   bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Parameters
  int N = 5.87 * mi.get_size();
  double rho = 0.02;
  double alpha = 1.0;
  int k = 10;
  int K = 1000;

  // Wrap the entire algorithm in a random restart (since p will eventually
  // converge, we don't expect to gain any benefit from running indefinitely).
  for (int iter=0; MaxCutHeuristic::Report(iter); ++iter) {
    // Alg 1 Step 1-3: Initialize probs of each index, best solution, and loop
    // vars
    std::vector<double> p(mi.get_size(), 0.5);
    double bestSolution = 0.0;
    int tprime = 0;  // Iterations since improvement
    int t = 0;  // Total iterations
    
    // Loop until termination criterion (or out of time)
    while (tprime < k && t < K) {
      // Alg 1 Step 4-5: Generate N random solutions using probabilities p,
      // keeping the best-performing rho*N. To do this in a memory-efficient way,
      // we will use a min heap, discarding solutions that are not among the best.
      int numKeep = std::max<int>(1, (int)(rho * N));
      std::vector<MaxCutSolution> X;
      for (int ct=0; ct < numKeep; ++ct) {
        X.push_back(MaxCutSolution::RandomSolution(mi, p, this));

        // Report each solution, as it could be a new best
        if (!Report(X[X.size()-1], iter)) {
          return;  // Out of time
        }
      }
      std::make_heap(X.begin(), X.end(), std::greater<MaxCutSolution>());
      for (int ct=0; ct < N-numKeep; ++ct) {
        MaxCutSolution newSolution = MaxCutSolution::RandomSolution(mi, p, this);
        bool added = false;
        if (newSolution.ImprovesOver(X.front())) {
          // Remove the worst solution from the heap and add this one
          std::pop_heap(X.begin(), X.end());
          X.pop_back();
          X.push_back(newSolution);
          std::push_heap(X.begin(), X.end());
          added = true;
        }
        // Report a solution if it was added (it could be a new best) and also
        // periodically (since this loop can take a while)
        if (added || ct % 10 == 0) {
          if (!Report(newSolution, iter)) {
            return;  // Out of time
          }
        }
      }

      // Alg 1 Steps 6-7: Update the p vector based on the best elements in X
      for (int j=0; j < mi.get_size(); ++j) {
	int ct = 0;
	for (int rank=0; rank < numKeep; ++rank) {
	  ct += (X[rank].get_assignments()[j] == 1);
	}
	p[j] = alpha * ct / numKeep + (1.0 - alpha) * p[j];
      }
      
      // Alg 1 Steps 8-10: Update loop vars
      double thisBest = X[0].get_weight();
      for (int j=1; j < X.size(); ++j) {
	thisBest = std::max<double>(thisBest, X[j].get_weight());
      }
      if (BaseSolution::ImprovesOver(thisBest, bestSolution)) {
	bestSolution = thisBest;
	tprime = 0;
      } else {
	++tprime;
      }
      ++t;
    }
  }
}

void Laguna2009HCE::LocallyOptimize(std::vector<FirstFixedMaxCutSolution>* X) {
  // Parameters
  double delta = 0.9;

  // Alg 2 Step 1: Initialize the state for the search; inL says if a node is
  // in the set L and distL is the sum of the hamming distances from each node
  // to all nodes in L.
  std::vector<bool> inL(X->size(), false);
  int numInL = 0;
  std::vector<int> distL(X->size(), 0);

  // First, identify the best solution in the population.
  int bestIdx = 0;
  double bestWeight = (*X)[0].get_weight();
  for (int idx=1; idx < X->size(); ++idx) {
    if ((*X)[idx].get_weight() > bestWeight) {
      bestIdx = idx;
      bestWeight = (*X)[idx].get_weight();
    }
  }
  
  // Alg 2 Step 2: Locally optimize the first solution and add it to L
  (*X)[bestIdx].AllShuffle1Swap();
  inL[bestIdx] = true;
  ++numInL;
  for (int i=0; i < X->size(); ++i) {
    if (!inL[i]) {
      distL[i] += (*X)[bestIdx].SymmetricDifference((*X)[i]);
    }
  }

  // Loop while |L| < delta * N
  while (numInL < delta * X->size()) {
    // Alg 2 Step 3: Find solution s furthest from L
    int bestDist = -1;
    int s = -1;
    for (int i=0; i < X->size(); ++i) {
      if (!inL[i] && distL[i] > bestDist) {
        s = i;
	bestDist = distL[i];
      }
    }

    // Alg 2 Step 4: Locally optimize s and add it to L
    (*X)[s].AllShuffle1Swap();
    inL[s] = true;
    ++numInL;
    if (numInL < delta * X->size()) {
      // If it's not the last iteration, increment distL values
      for (int i=0; i < X->size(); ++i) {
	if (!inL[i]) {
	  distL[i] += (*X)[s].SymmetricDifference((*X)[i]);
	}
      }
    }
  }
}

Laguna2009HCE::Laguna2009HCE(const MaxCutInstance& mi, double runtime_limit,
			     bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Parameters
  int N = std::max<int>(1, 0.031 * mi.get_size());
  double alpha = 0.9;
  int K = 100;
  
  // Wrap entire algorithm with specified K parameter into a random restart;
  // since the p vector eventually converges running it indefinitely is not
  // likely to improve solutions.
  for (int iter=0; MaxCutHeuristic::Report(iter); ++iter) {
    // Alg 3 Step 1: Initialize vertex probabilities to all 0.5
    std::vector<double> p(mi.get_size(), 0.5);
    for (int t=0; t < K; ++t) {
      // Alg 3 Step 4: Generate N random solutions using probabilities p
      std::vector<FirstFixedMaxCutSolution> X;
      while (X.size() < N) {
	X.push_back(FirstFixedMaxCutSolution::RandomSolution(mi, p, this, 1));
      }
      
      // Alg 3 Steps 5-6: Apply local optimization procedure from Fig. 2
      LocallyOptimize(&X);
      
      // Alg 3 Steps 7-8: Update p vector based on all elements of X
      for (int j=0; j < mi.get_size(); ++j) {
	int ct = 0;
	for (int idx=0; idx < N; ++idx) {
	  ct += (X[idx].get_assignments()[j] == 1);
	}
	p[j] = alpha * ct / N + (1.0 - alpha) * p[j];
      }
      
      // Check termination criterion for best solution in X
      double bestWt = X[0].get_weight();
      int bestIdx = 0;
      for (int idx=1; idx < X.size(); ++idx) {
	if (X[idx].get_weight() > bestWt) {
          bestWt = X[idx].get_weight();
	  bestIdx = idx;
	}
      }
      if (!Report(X[bestIdx], iter)) {
	return;  // Out of time
      }
    }
  }
}
