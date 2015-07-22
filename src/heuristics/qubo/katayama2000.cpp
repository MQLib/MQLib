#include <algorithm>
#include <iostream>
#include "heuristics/qubo/katayama2000.h"
#include "util/random.h"

Katayama2000QUBOSolution::Katayama2000QUBOSolution(const Katayama2000QUBOSolution &x1,
						   const Katayama2000QUBOSolution &x2) :
  QUBOSolution(x1) {
  // Parameters
  int mutateDist = N_/10;

  // If both solutions have same value, use that. Otherwise, randomly assign.
  int dParents = 0;  // Hamming distance between x1 and x2
  int d1 = 0;  // Distance of new solution from x1
  int d2 = 0;  // Distance of new solution from x2
  for (int i=0; i < N_; ++i) {
    if (x1.assignments_[i] != x2.assignments_[i]) {
      ++dParents;
      if (Random::RandDouble() < 0.5) {
	// Keep the x1 value
	++d2;
      } else {
	// Switch to the x2 value
	++d1;
	UpdateCutValues(i);
      }
    }
  }

  // If the hamming distance between the parents is less than mutateDist, flip
  // randomly selected indices where the two solutions are the same until the
  // dist from this solution to one of the parents is mutateDist.
  if (dParents < mutateDist) {
    int numFlip = std::min<int>(mutateDist-d1, mutateDist-d2);
    std::vector<int> identical;
    for (int i=0; i < N_; ++i) {
      if (x1.assignments_[i] == x2.assignments_[i]) {
	identical.push_back(i);
      }
    }
    std::random_shuffle(identical.begin(), identical.end());
    for (int i=0; i < numFlip; ++i) {
      UpdateCutValues(identical[i]);
    }
  }
}

void Katayama2000QUBOSolution::VariantKOpt() {
  // best is the best solution found in the LS. At the end of each iteration,
  // we will copy over this best solution to be the current solution.
  Katayama2000QUBOSolution best(*this);

  // Loop until no improvement during loop
  bool improved = true;
  while (improved) {
    // Initialize iteration variables (in terminology of Figure 1, improved
    // becomes true when Gmax > 0).
    improved = false;
    std::vector<bool> inC(N_, true);  // Is a variable in set C (not flipped)?
    int numInC = N_;

    // Loop until no variables are in C
    while (numInC > 0) {
      // Fig 1 Step 1.2.1: Generate random permutation of variables
      std::vector<int> RP(N_, 0);
      for (int i=0; i < N_; ++i) {
	RP[i] = i;
      }
      std::random_shuffle(RP.begin(), RP.end());

      // Fig 1 Step 1.2.2: Search all variables (regardless of whether they're
      // in C) in RP order, and flip a variable if doing so will improve best.
      for (int i=0; i < N_; ++i) {
	if (ImprovesOverAfterMove(best, RP[i])) {
	  improved = true;
	  UpdateCutValues(RP[i]);
	  if (inC[RP[i]]) {
	    inC[RP[i]] = false;
	    --numInC;
	  }
	  best = *this;
	}
      }

      // Break out of loop if we've now exhausted C
      if (numInC == 0) {
	break;
      }

      // Fig 1 Steps 1.2.3: Find index j (j \in C) with best diff_weights_[j]
      int j = -1;
      double bestGain = -std::numeric_limits<double>::max();
      for (int i=0; i < N_; ++i) {
	if (inC[i] && diff_weights_[i] > bestGain) {
	  j = i;
	  bestGain = diff_weights_[i];
	}
      }

      // Fig 1 Steps 1.2.4, 1.2.6-1.2.7: Flip index j
      UpdateCutValues(j);

      // Fig 1 Step 1.2.5: If new solution improves on best, copy it over
      if (ImprovesOver(best)) {
	improved = true;
	best = *this;
      }

      // Fig 1 Step 1.2.8: Remove j from set C
      inC[j] = false;
      --numInC;
    }

    // Fig 1 Step 1.3: Copy over best to be the current solution
    operator=(best);
  }
}

void Katayama2000QUBOSolution::Mutate() {
  // Parameters
  int numFlip = N_ / 2;

  // Randomly select the numFlip indices to flip
  std::vector<int> toflip(N_, 0);
  for (int i=0; i < N_; ++i) {
    toflip[i] = i;
  }
  std::random_shuffle(toflip.begin(), toflip.end());

  // Flip selected bits
  for (int i=0; i < numFlip; ++i) {
    UpdateCutValues(toflip[i]);
  }
}

Katayama2000Elite::Katayama2000Elite(const QUBOInstance& qi, int PS,
				     QUBOHeuristic *heuristic) :
  PS_(PS),
  heuristic_(heuristic),
  stepsSinceImprovement_(0) {
  // Initialize a random population of size PS, and run local search on each
  for (int ct=0; ct < PS; ++ct) {
    P_.push_back(QUBOSolution::RandomSolution(qi, heuristic));
    P_[ct].VariantKOpt();

    // This is a pretty time-consuming initialization, so we'll report each new
    // solution.
    if (!heuristic->Report(P_[ct])) {
      return;  // Out of time
    }
  }

  // Sort P_ by objective value
  std::sort(P_.begin(), P_.end(), std::greater<Katayama2000QUBOSolution>());
}

void Katayama2000Elite::Update(const std::vector<Katayama2000QUBOSolution>& x) {
  // Obtain old best objective (used to determine if we've improved)
  double oldBest = P_[0].get_weight();

  // Build a combined list of solutions and sort by objective
  std::vector<Katayama2000QUBOSolution> combined(P_);
  combined.insert(combined.end(), x.begin(), x.end());
  std::sort(combined.begin(), combined.end(),
  	    std::greater<Katayama2000QUBOSolution>());

  // Add the best non-duplicated solutions to P_.
  P_.clear();
  for (int i=0; i < combined.size(); ++i) {
    // Check if it's duplicated
    bool match = false;
    for (int j=0; j < P_.size(); ++j) {
      if (combined[i] == P_[j]) {
	match = true;
	break;
      }
    }

    if (!match) {
      P_.push_back(combined[i]);
    }
    
    if (P_.size() == PS_) {
      break;  // We have selected the maximum number of elements
    }
  }

  // Update stepsSinceImprovement_
  if (P_[0].ImprovesOver(oldBest)) {
    stepsSinceImprovement_ = 0;
  } else {
    ++stepsSinceImprovement_;
  }
}

void Katayama2000Elite::Diversify() {
  // Parameters
  int diversitySteps = 30;  // Steps without improvement before diversification
  double avgDistThreshold = 30.0;  // Min average pairwise dist before diversify

  bool needDiversify = false;
  if (stepsSinceImprovement_ >= diversitySteps) {
    needDiversify = true;
  } else {
    // Check pairwise hamming distances
    int hammingSum = 0;
    for (int i=0; i < P_.size(); ++i) {
      for (int j = i+1; j < P_.size(); ++j) {
	hammingSum += P_[i].SymmetricDifference(P_[j]);
      }
    }
    if (hammingSum < avgDistThreshold * P_.size() * (P_.size() - 1) / 2.0) {
      needDiversify = true;
    }
  }

  // If we need diversification, randomly select n/2 bits and flip them,
  // followed by local search. Do this for all except best solution. Also reset
  // stepsSinceImprovement_ in this case.
  if (needDiversify) {
    // Reset stepsSinceImprovement_
    stepsSinceImprovement_ = 0;

    // Mutate and local search on each solution
    for (int i=1; i < P_.size(); ++i) {
      P_[i].Mutate();
      P_[i].VariantKOpt();

      // This loop could take a while, so report solutions in the loop
      if (!heuristic_->Report(P_[i])) {
        return;
      }
    }

    // Re-sort population by objective
    std::sort(P_.begin(), P_.end(), std::greater<Katayama2000QUBOSolution>());
  }
}

Katayama2000::Katayama2000(const QUBOInstance& qi, double runtime_limit,
			   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  int PS = 20;

  // Fig 2 Steps 1-2: Initialize a random population of size PS, and run
  // local search on each
  Katayama2000Elite P(qi, PS, this);

  // Exit if we're already out of time
  if (!Report(P.get_best())) {
    return;
  }

  // Fig 2 Step 3: Loop until termination criterion (runtime limit)
  while (1) {
    // Fig 2 Step 3.2: Build PS/2 offspring
    std::vector<Katayama2000QUBOSolution> offspring;
    for (int ct=0; ct < PS/2; ++ct) {
      // Fig 2 Step 3.2.1: Select random parents. While we specify PS to be the
      // size of the population, sometimes it falls below this size (in
      // problems where local search from random solutions almost always
      // yield identical solutions). Therefore, we need to check the size of
      // P instead of assuming it is PS.
      int currPS = P.get_size();
      int a, b;  // Parent indices
      if (currPS == 1) {
	a = 0;
	b = 0;
      } else {
	do {
	  a = Random::RandInt(0, currPS-1);
	  b = Random::RandInt(0, currPS-1);
	} while (a == b);
      }

      // Fig 2 Step 3.2.2-3.2.3: Obtain offspring via crossover, performing
      // mutation if necessary.
      offspring.push_back(Katayama2000QUBOSolution::Crossover(P.get_P()[a],
							      P.get_P()[b]));

      // Fig 2 Step 3.2.4: Perform local search on new offspring
      offspring[ct].VariantKOpt();

      // The variant k-opt is expensive, so we'll report each child
      if (!Report(P.get_best())) {
        break;
      }
    }

    // Fig 2 Step 3.3: Update the elite set of solutions
    P.Update(offspring);

    // Check termination criterion with the best elite solution
    if (!Report(P.get_best())) {
      break;
    }

    // Fig 2 Step 3.4: Perform diversification strategy if needed
    P.Diversify();

    // Check termination criterion again, as Diversify can take a long time
    if (!Report(P.get_best())) {
      break;
    }
  }
}
