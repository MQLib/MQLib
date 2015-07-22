#include <algorithm>
#include <iostream>
#include "heuristics/qubo/merz2004.h"
#include "util/random.h"

// Crossover -- we'll start with x1 and just flip the variables that change
Merz2004Solution::Merz2004Solution(const Merz2004Solution& x1,
				   const Merz2004Solution& x2) :
  QUBOSolution(x1) {
  // Fig 7 steps 2-7 build sets of indices different and the same between x1, x2
  std::vector<int> nCB;  // Indices that differ between x1 and x2
  std::vector<int> CB;  // Indices that are common between x1 and x2
  int num = x1.SymmetricDifference(x2, &nCB, &CB);

  // Fig 7 step 8: Repeat loop "num" times
  for (int t=0; t < num; ++t) {
    // Fig 7 steps 9-12: Flip a random non-common bit with positive diff_weights_
    std::vector<int> positive;  // Positions in nCB of positive diff_weights_
    for (int i=0; i < nCB.size(); ++i) {
      if (ImprovingMove(nCB[i])) {
	positive.push_back(i);
      }
    }

    if (positive.size() > 0) {
      // Fig 7 Step 9: Randomly select the variable
      int idx = positive[Random::RandInt(0, positive.size()-1)];  // pos in nCB

      // Fig 7 Steps 10-11: Flip the variable and update gains
      UpdateCutValues(nCB[idx]);

      // Fig 7 Step 12: Remove variable from nCB vector
      nCB[idx] = nCB[nCB.size()-1];
      nCB.resize(nCB.size()-1);
    }

    // Fig 7 Steps 14-18: If CB is non-empty, flip the best variable in it
    if (CB.size() > 0) {
      // Fig 7 Step 14: Find the element of CB with the best gain
      int bestIdx = 0;
      double bestGain = diff_weights_[CB[0]];
      for (int i=1; i < CB.size(); ++i) {
	if (diff_weights_[CB[i]] > bestGain) {
	  bestIdx = i;
	  bestGain = diff_weights_[CB[i]];
	}
      }

      // Fig 7 Steps 15-16: Flip CB[bestIdx] and update gains
      UpdateCutValues(CB[bestIdx]);

      // Fig 7 Step 17: Remove flipped variable from CB
      CB[bestIdx] = CB[CB.size()-1];
      CB.resize(CB.size()-1);
    }
  }
}

void Merz2004Solution::RandomizedKOpt() {
  // Parameters
  int m = 50;  // Number of inner loop iterations without improvement before exit

  // best is the best solution found in the LS. At the end of each iteration,
  // we will copy over this best solution to be the current solution.
  Merz2004Solution best(*this);

  // Loop until no improvement during loop
  bool improved = true;
  while (improved) {
    // Initialize iteration variables (in terminology of Figure 1, improved
    // becomes true when Gmax > 0).
    improved = false;
    std::vector<bool> inC(N_, true);  // Is a variable in set C (not flipped)?
    int numInC = N_;

    // Loop until no variables are in C
    int stepsSinceImprovement = 0;
    while (stepsSinceImprovement <= m && numInC > 0) {
      // Increment steps since improvement
      ++stepsSinceImprovement;

      // Fig 1 Step 1.2.1: Generate random permutation of variables
      std::vector<int> RP;
      for (int i=0; i < N_; ++i) {
	RP.push_back(i);
      }
      std::random_shuffle(RP.begin(), RP.end());

      // Fig 1 Step 1.2.2: Search all variables (regardless of whether they're
      // in C) in RP order, and flip a variable if doing so will improve best.
      for (int i=0; i < N_; ++i) {
	if (ImprovesOverAfterMove(best, RP[i])) {
	  improved = true;
	  stepsSinceImprovement = 0;
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
	stepsSinceImprovement = 0;
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

void Merz2004Solution::Mutate() {
  // Parameters
  int numFlip = N_ / 3;
  
  // Randomly select the numFlip indices to flip
  std::vector<int> toflip;
  for (int i=0; i < N_; ++i) {
    toflip.push_back(i);
  }
  std::random_shuffle(toflip.begin(), toflip.end());

  // Flip selected bits
  for (int i=0; i < numFlip; ++i) {
    UpdateCutValues(toflip[i]);
  }
}

Merz2004Elite::Merz2004Elite(const QUBOInstance& qi, int PS,
			     QUBOHeuristic *heuristic) :
  PS_(PS) {
  // Initialize a random population of size pS, and run local search on each. Sort
  // by objective value.
  std::vector<Merz2004Solution> initial;
  for (int ct=0; ct < PS; ++ct) {
    initial.push_back(Merz2004Solution::RandomizedGreedy(qi, heuristic));
    initial[ct].RandomizedKOpt();
    if (!heuristic->Report(initial[ct])) {
      break;  // Stop generating new initial solutions if runtime limit hit
    }
  }

  // Add the non-duplicated solution to P_ (note that we may end up with fewer than
  // PS_ total solutions).
  SelectNonDuplicated(&initial);
}

void Merz2004Elite::SelectNonDuplicated(std::vector<Merz2004Solution>* x) {
  // Sort x by objective value
  std::sort(x->begin(), x->end(), std::greater<Merz2004Solution>());

  // Fill P_ with the non-duplicated values from x, but select no more than PS_
  // total values (note that there may be less than PS_ final values if there are
  // duplicates).
  P_.clear();
  for (int i=0; i < x->size(); ++i) {
    // Check if it's duplicated
    bool match = false;
    for (int j=0; j < P_.size(); ++j) {
      if ((*x)[i] == P_[j]) {
	match = true;
	break;
      }
    }

    if (!match) {
      P_.push_back((*x)[i]);
    }
    
    if (P_.size() == PS_) {
      break;  // We have selected all the elements
    }
  }
}

void Merz2004Elite::Update(const std::vector<Merz2004Solution>& x) {
  // Obtain old best objective (used to determine if we've improved)
  double oldBest = P_[0].get_weight();

  // Build a combined list of solutions and sort by objective
  std::vector<Merz2004Solution> combined(P_);
  combined.insert(combined.end(), x.begin(), x.end());

  // Add the best non-duplicated solutions to P_.
  SelectNonDuplicated(&combined);

  // Update stepsSinceImprovement_
  if (P_[0].ImprovesOver(oldBest)) {
    stepsSinceImprovement_ = 0;
  } else {
    ++stepsSinceImprovement_;
  }
}

void Merz2004Elite::Diversify() {
  // Parameters
  int diversitySteps = 30;  // Steps without improvement before diversification

  // If we need diversification, mutate then local search all solutions except
  // the best. Also reset stepsSinceImprovement_.
  if (stepsSinceImprovement_ >= diversitySteps) {
    // Reset stepsSinceImprovement_
    stepsSinceImprovement_ = 0;

    // Mutate and local search on each solution except the best (position 0)
    for (int i=1; i < P_.size(); ++i) {
      P_[i].Mutate();
      P_[i].RandomizedKOpt();
    }

    // Re-sort population by objective
    std::sort(P_.begin(), P_.end(), std::greater<Merz2004Solution>());
  }
}

Merz2004::Merz2004(const QUBOInstance& qi, double runtime_limit,
		   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  int PS = 40;
  double crossoverRate = 0.5;

  // Fig 5 Steps 1-2: Initialize population and perform local search on each
  Merz2004Elite P(qi, PS, this);
  if (!Report(P.get_best())) {
    return;
  }

  // Fig 5 Step 3: Loop until termination criterion (we will break when out of
  //               time)
  while (true) {
    // Fig 5 Steps 4-9: Generate new offspring in a set Pc
    std::vector<Merz2004Solution> Pc;
    for (int idx=0; idx < crossoverRate * PS; ++idx) {
      // Fig 5 Step 5: Select two parents at random from P
      int popSize = P.get_P().size();
      int a, b;
      if (popSize == 1) {
        a = 0;
        b = 0;
      } else {
        do {
          a = Random::RandInt(0, popSize-1);
          b = Random::RandInt(0, popSize-1);
        } while (a == b);
      }

      // Fig 5 Step 6: Generate new offspring via crossover
      Pc.push_back(Merz2004Solution::Crossover(P.get_P()[a], P.get_P()[b]));

      // Fig 5 Step 7: Perform local search on the new solution
      Pc[idx].RandomizedKOpt();
    }

    // Fig 5 Step 10: Update the elite set of solutions
    P.Update(Pc);

    // Check termination criterion with best elite solution
    if (!Report(P.get_best())) {
      break;
    }

    // Fig 5 Steps 11-13: Perform diversification strategy if needed
    P.Diversify();
  }
}
