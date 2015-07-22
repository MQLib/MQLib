#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/qubo/pardalos2008.h"
#include "util/random.h"

Pardalos2008QUBOSolution::Pardalos2008QUBOSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

Pardalos2008QUBOSolution::Pardalos2008QUBOSolution
(const Pardalos2008QUBOSolution x, const std::vector<double>& probs, int k) :
  QUBOSolution(x) {
  // Parameters
  int distmax = (k <= 5) ? N_/2 : N_/5;

  // Alg 3 Step 1: Initialize dist, the number of values flipped
  int dist = 0;

  // Alg 3 Steps 2, 12: Loop through variables in order
  for (int j=0; j < N_; ++j) {
    // Alg 3 Steps 3-11: Probabilistically flip nodes
    if (assignments_[j]) {
      if (probs[j] <= Random::RandDouble()) {
	UpdateCutValues(j);
	++dist;
      }
    } else {
      if (probs[j] >= Random::RandDouble()) {
	UpdateCutValues(j);
	++dist;
      }
    }

    // Alg 3 Steps 13-15: If we have already flipped distmax, stop.
    if (dist == distmax) {
      break;
    }
  }
}

void Pardalos2008QUBOSolution::TabuSearch
(std::vector<Pardalos2008QUBOSolution>* R) {
  // Parameters
  int nbad = N_;  // Tabu search stopping parameter
  int tabu = 21;  // Tabu tenure

  // Alg 2 Step 1: Initialize parameters
  Pardalos2008QUBOSolution xbest = *this;
  std::vector<int> M;
  for (int idx=0; idx < N_; ++idx) {
    M.push_back(idx);
  }
  int step = 0;
  bool impr = true;
  R->clear();

  // Alg 2 Step 2: Loop until no improvement
  while (impr) {
    // Alg 2 Step 3: Initialize vars specific to each iteration
    impr = false;
    double deltaBest = 0.0;
    double deltaCur = 0.0;
    int stepImpr = step;
    std::vector<int> last_used(N_, -tabu - 1);
    
    // Alg 2 Steps 4, 39: Loop until step-stepImpr > nbad (Step 39 says the
    // opposite, but it seems Fig 2 used "until" instead of "while").
    while (step - stepImpr < nbad) {
      // Alg 2 Steps 5, 21: Add all improving moves until none exist, making
      // them in random order. Again, Step 21 says the opposite, but that would
      // cause an infinite loop for any locally optimal x. As a result, it seems
      // "until" there should have instead been "while"
      double delta;
      do {
	delta = 0.0;

	// Alg 2 Steps 7-9: If the solution is locally optimal in the 1-flip
	// neighborhood (all 1-flips are non-improving), save it in R.
	bool allNonImprove = true;
	for (int idx=0; idx < N_; ++idx) {
	  if (ImprovingMove(idx)) {
	    allNonImprove = false;
	    break;
	  }
	}
	if (allNonImprove) {
	  R->push_back(*this);
	}

	// Alg 2 Step 10: Randomly permute M
	std::random_shuffle(M.begin(), M.end());

	// Alg 2 Steps 11-20: Go through the variables in the order of M,
	// flipping if they're non-decreasing 1-moves and either allowed by the
	// tabu list or will improve on xbest.
	for (int k=0; k < N_; ++k) {
	  int j = M[k];
	  if (step - last_used[j] > tabu || ImprovesOverAfterMove(xbest, j)) {
	    if (ImprovingMove(j)) {
	      delta += diff_weights_[j];
	      last_used[j] = step;
	      UpdateCutValues(j);
	      ++step;
	    }
	  }
	}
	deltaCur += delta;
      } while (delta > 0.0);

      // Alg 2 Steps 22-27: If the current solution improves on or matches
      // xbest, replace it. It it's improving, note this in stepImpr and impr.
      if (!BaseSolution::ImprovesOver(deltaBest, deltaCur)) {
	if (BaseSolution::ImprovesOver(deltaCur, deltaBest)) {
	  deltaBest = deltaCur;
	  stepImpr = step;
	  impr = true;
	}
	xbest = *this;
      }

      // Alg 2 Steps 28-36: Identify the best one-flip not blocked by the
      // tabu tenure (regardless of whether it's improving or not).
      delta = -std::numeric_limits<double>::max();
      int indBest = -1;
      for (int k=0; k < N_; ++k) {
	int j = M[k];
	if (step - last_used[j] > tabu || ImprovesOverAfterMove(xbest, j)) {
	  delta = diff_weights_[j];
	  indBest = j;
	}
      }

      // Alg 2 Steps 37-38: Perform the 1-flip we just identified (even if it
      // makes us worse)
      if (indBest >= 0) {
        deltaCur += delta;
        last_used[indBest] = step;
        UpdateCutValues(indBest);
      }
      ++step;
    }

    // Alg 2 Step 40: Set this solution to be the best ever identified
    operator=(xbest);
  }
}

Pardalos2008Elite::Pardalos2008Elite(int Esize) :
  Esize_(Esize) {}

const Pardalos2008QUBOSolution& Pardalos2008Elite::Best() const {
  if (Elite_.size() == 0) {
    std::cout << "Called Best() without any elite solutions" << std::endl;
    exit(0);
  }
  int best_idx = 0;
  double best_weight = Elite_[0].get_weight();
  for (int idx=1; idx < Elite_.size(); ++idx) {
    if (Elite_[idx].get_weight() > best_weight) {
      best_idx = idx;
      best_weight = Elite_[idx].get_weight();
    }
  }
  return Elite_[best_idx];
}


void Pardalos2008Elite::AddSolution(const Pardalos2008QUBOSolution& x) {
  std::vector<Pardalos2008QUBOSolution> toAdd;
  toAdd.push_back(x);
  AddSolutions(toAdd);
}

void Pardalos2008Elite::AddSolutions(const std::vector<
				     Pardalos2008QUBOSolution>& slns) {
  // Check if each solution can be added to the list of elite solutions
  for (int idx=0; idx < slns.size(); ++idx) {
    if (Elite_.size() < Esize_) {
      // Elite_ is not full; just add the new solution and re-make heap
      Elite_.push_back(slns[idx]);
      std::push_heap(Elite_.begin(), Elite_.end(),
		     std::greater<Pardalos2008QUBOSolution>());
    } else if (slns[idx].ImprovesOver(Elite_[0])) {
      // Elite_ is full but new solution is better than the worst in Elite_,
      // so replace it
      std::pop_heap(Elite_.begin(), Elite_.end(),
		    std::greater<Pardalos2008QUBOSolution>());
      Elite_.pop_back();
      Elite_.push_back(slns[idx]);
      std::push_heap(Elite_.begin(), Elite_.end(),
		     std::greater<Pardalos2008QUBOSolution>());
    }
  }
}

void Pardalos2008Elite::LimitByBests(const std::vector<
				     Pardalos2008QUBOSolution>& bests) {
  // Parameters
  int dp = 200;  // Prohibition parameter

  // If there are fewer than dp variables, just throw away all elite solutions.
  if (Elite_.size() > 0 && Elite_[0].get_assignments().size() < dp &&
      bests.size() > 0) {
    Elite_.clear();
    return;
  }

  // Remove any elite solution that is within dp of any of the passed best
  // solutions.
  std::vector<Pardalos2008QUBOSolution> NewElite;
  for (int i=0; i < Elite_.size(); ++i) {
    bool tooClose = false;
    for (int j=0; j < bests.size(); ++j) {
      if (Elite_[i].SymmetricDifference(bests[j]) <= dp) {
	tooClose = true;
	break;
      }
    }
    if (!tooClose) {
      NewElite.push_back(Elite_[i]);
    }
  }
  std::make_heap(NewElite.begin(), NewElite.end(),
		 std::greater<Pardalos2008QUBOSolution>());
  Elite_ = NewElite;
}

Pardalos2008Probs::Pardalos2008Probs(const std::vector<
				     Pardalos2008QUBOSolution>& slns, int K,
				     const std::vector<double>& mu) :
  K_(K),
  mu_(mu),
  N_(slns[0].get_assignments().size()),
  numerator0_((K_+1)*N_, 0.0),
  numerator1_((K_+1)*N_, 0.0),
  denominator0_((K_+1)*N_, 0.0),
  denominator1_((K_+1)*N_, 0.0),
  freq1_(N_, 0),
  freq_(0) {
  AddSolutions(slns);
}

void Pardalos2008Probs::AddSolutions(const std::vector<
				     Pardalos2008QUBOSolution>& slns) {
  // Update frequency information
  freq_ += slns.size();
  for (int idx=0; idx < slns.size(); ++idx) {
    for (int j=0; j < N_; ++j) {
      freq1_[j] += (slns[idx].get_assignments()[j] == 1);
    }
  }

  // Update numerator and denominator information
  for (int idx=0; idx < slns.size(); ++idx) {
    for (int k=0; k <= K_; ++k) {
      double expWeight = exp(-mu_[k] * slns[idx].get_weight());
      for (int j=0; j < N_; ++j) {
	if (slns[idx].get_assignments()[j]) {
	  numerator1_[k*N_ + j] += slns[idx].get_weight() * expWeight;
	  denominator1_[k*N_ + j] += expWeight;
	} else {
	  numerator0_[k*N_ + j] += slns[idx].get_weight() * expWeight;
	  denominator0_[k*N_ + j] += expWeight;
	}
      }
    }
  }
}

void Pardalos2008Probs::GetProbs(int k, std::vector<double>* probs) const {
  probs->clear();
  for (int j=0; j < N_; ++j) {
    if (freq1_[j] == 0) {
      probs->push_back(0.0);
    } else if (freq1_[j] == freq_) {
      probs->push_back(1.0);
    } else {
      // We have seen both a 1 and 0 in position j; compute the prob
      double pj0 = ((double)freq1_[j]) / freq_;  // Prop of position j as 1
      double expval = 0.0;
      for (int i=0; i < k; ++i) {
	expval -= 0.5 * (mu_[i+1] - mu_[i]) *
	  (numerator0_[i*N_ + j] / denominator0_[i*N_ + j] +
	   numerator0_[(i+1)*N_ + j] / denominator0_[(i+1)*N_ + j] -
	   numerator1_[i*N_ + j] / denominator1_[i*N_ + j] -
	   numerator1_[(i+1)*N_ + j] / denominator1_[(i+1)*N_ + j]);
      }
      probs->push_back(1.0 / (1.0 + (1.0 - pj0) / pj0 * exp(expval)));
    }
  }
}


Pardalos2008::Pardalos2008(const QUBOInstance& qi, double runtime_limit,
			   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  int K = 25;  // Number of temperature stages
  std::vector<double> mu(K+1, 0.0);  // Temperature schedule
  mu[0] = 0.0;
  mu[1] = 0.0000001;
  for (int k=2; k <= K; ++k) {
    mu[k] = mu[k-1] * ((log(10.0) - log(mu[1])) / 48.0);
  }
  int ngen = 27;  // Number of solutions generated
  int maxnfail = 1;  // Max failures before exit of middle loop
  bool restartCriterion = false;  // True only for maximum independent set pbms
  int Esize = qi.get_size() / 2;  // Max number of elite solutions

  // Initialize variables
  Pardalos2008Elite Elite(Esize);  // Set of elite solutions
  Pardalos2008QUBOSolution xbest = QUBOSolution::RandomSolution(qi, this);
  std::vector<Pardalos2008QUBOSolution> bests;  // Set of iter.-end best slns

  // Alg 1 step 2: Loop until termination criterion met
  while (1) {
    // Alg 1 steps 3-8: Stilde is initialized to be equal to Elite in all cases;
    // if Elite starts as empty then we'll need to initialize it. This also
    // causes P (set of neighborhoods around best solutions) to be emptied.
    if (Elite.size() == 0) {
      Elite.AddSolution(QUBOSolution::RandomSolution(qi, this));
      bests.clear();
    }

    // xmax is initialized to the best elite solution at the start of each
    // iteration (Alg 1, steps 3, 35)
    Pardalos2008QUBOSolution xmax = Elite.Best();

    // Alg 1, Step 9: Repeat main loop until maxnfail iterations w/o improvement
    int nfail = 0;
    while (nfail <= maxnfail) {
      // Alg 1, Step 10: Remember objective of xmax at beginning of iteration
      double fxold = xmax.get_weight();

      // Stilde is actually incorporated into a Pardalos2008Probs data structure
      // that maintains the numerator and denominator of E_{kj}^u. For each k,
      // j, and u. It is initialized each iteration from Elite (Alg. 1 steps 6,
      // 27, and 34)
      Pardalos2008Probs probs(Elite.get_solutions(), K, mu);

      // Alg 1, Step 11: Loop through temperature cycles
      for (int k=0; k <= K; ++k) {
	// Alg 1, Step 12: Calculate generation probabilities
	std::vector<double> ptilde;
	probs.GetProbs(k, &ptilde);

	// Alg 1, Step 13: Repeat generation steps ngen+1 times
	for (int g=0; g <= ngen; ++g) {
	  // Alg 1, Step 14: Generate a new solution from xmax and ptilde
	  Pardalos2008QUBOSolution x =
	    Pardalos2008QUBOSolution::GenerateSolution(xmax, ptilde, k);

	  // Alg 1, Step 15: Get a set of locally optimal solutions via tabu
	  // search
	  std::vector<Pardalos2008QUBOSolution> R;
	  x.TabuSearch(&R);

	  // Alg 1, Steps 16-21: Update probs, Elite, xmax, and xbest based on
	  // the set R.
	  probs.AddSolutions(R);
	  Elite.AddSolutions(R);
	  for (int idx=0; idx < R.size(); ++idx) {
	    if (R[idx].ImprovesOver(xmax)) {
	      xmax = R[idx];
	    }
	    if (R[idx].ImprovesOver(xbest)) {
	      xbest = R[idx];
	    }
	  }

          // Let's check our termination criterion here
          if (!Report(xbest)) {
            return;  // Runtime limit reached
          }
	}
      }
      
      // Alg 1, Steps 24-26: Check if xmax improved this iteration
      if (xmax.ImprovesOver(fxold)) {
	nfail = 0;
      } else {
	++nfail;
      }
    }

    // Alg 1 Steps 29-33: Clear Elite if restartCriterion is true. Otherwise,
    // clear it based on recent best solutions.
    if (restartCriterion) {
      Elite.clear();
    } else {
      bests.push_back(xbest);
      Elite.LimitByBests(bests);
    }
  }
}
