#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "heuristics/qubo/glover2010.h"
#include "util/random.h"

Glover2010QUBOSolution::Glover2010QUBOSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

Glover2010QUBOSolution::Glover2010QUBOSolution(const Glover2010QUBOSolution& x,
					       const std::vector<int>&
					       EliteFreq, int r,
					       const std::vector<int>&
					       FlipFreq)
  : QUBOSolution(x) {
  // Parameters
  double beta = 0.3;  // Frequency-related weight in perturbation scoring
  double lambda = 1.2;  // Perturbation selection importance factor
  int gamma = N_ / 4;  // Perturbation strength

  // Step 1: Compute Score(x_i) from formula (3) for each variable. Create a
  // list of (score, index) pairs so we can sort and identify the indices.
  std::vector<std::pair<double, int> > Score;
  double max_Freq = (double)*std::max_element(FlipFreq.begin(), FlipFreq.end());

  for (int i=0; i < N_; ++i) {
    double score = ((double)(EliteFreq[i] * (r-EliteFreq[i]))) / r / r +
      beta * (1.0 - FlipFreq[i] / max_Freq);
    Score.push_back(std::pair<double, int>(score, i));
  }

  // Step 2: Sort the scores and compute the weight P for each variable (in the
  // paper these are normalized, but actually we'll use roulette wheel selection
  // so no need to normalize.
  std::sort(Score.begin(), Score.end(),
	    std::greater<std::pair<double, int> >());

  std::vector<double> P(N_);
  for (int i=0; i < N_; ++i) {
    P[Score[i].second] = pow(i + 1.0, -lambda);
  }

  // Step 3: Use multi-element roulette wheel selection with weights P to
  // select the gamma elements to flip. Section 2.3.2 seems to suggest exactly
  // gamma values are selected to flip, though the exact mechanism is not
  // worded precisely.
  std::vector<int> flips;
  Random::MultiRouletteWheel(P, gamma, &flips);

  // Step 4: Flip each selected index
  for (int idx=0; idx < flips.size(); ++idx) {
    UpdateCutValues(flips[idx]);
  }
}

void Glover2010QUBOSolution::TabuSearch(std::vector<int>* FlipFreq) {
  // Parameters
  int c = N_ / 100;  // Tabu tenure
  int randTenure = 10;  // Upper bound on random component of tabu tenure
  int alpha = 20 * N_;  // Improvement cutoff of TS

  // NonTabuTimes is the vector of the first iteration when each variable will
  // be non-tabu.
  std::vector<int> NonTabuTimes(N_, 0);

  // lastImprovement: step when we last improved
  int lastImprovement = 0;

  // Initialize FlipFreq, which tracks number of times we flip each var
  FlipFreq->clear();
  FlipFreq->resize(N_, 1);

  // We will store the best solution encountered during the search
  Glover2010QUBOSolution best(*this);

  int step = 0;  // Iteration count
  while (step - lastImprovement < alpha) {
    // Find the best non-tabu move (allow tabu moves if they will improve on the
    // best solution so far.
    double bestMove = -std::numeric_limits<double>::max();
    int bestIdx = -1;
    for (int i=0; i < N_; ++i) {
      if (NonTabuTimes[i] <= step || ImprovesOverAfterMove(best, i)) {
	if (diff_weights_[i] > bestMove) {
	  bestMove = diff_weights_[i];
	  bestIdx = i;
	}
      }
    }

    // Perform the selected move, and increment the associated tabu tenure
    if (bestIdx >= 0) {
      UpdateCutValues(bestIdx);
      NonTabuTimes[bestIdx] = step + c + Random::RandInt(1, randTenure);
      ++((*FlipFreq)[bestIdx]);
    }

    // If we've improved the best solution, copy it and update lastImprovement
    if (ImprovesOver(best)) {
      best = *this;
      lastImprovement = step;
    }

    // Increment step count
    ++step;
  }

  // Copy over the best solution in the tabu search as the final one
  operator=(best);
}

Glover2010Elite::Glover2010Elite(const QUBOInstance &qi, int R) :
  R_(R),
  N_(qi.get_size()),
  EliteFreq_(N_, 0) {}

void Glover2010Elite::AddSolution(const Glover2010QUBOSolution& x) {
  // First, we'll do the easy check: reject the new solution if it's not better
  // than the worst elite solution and we are at maximum capacity for the
  // elite set.

  if (EliteSol_.size() == R_ && !x.ImprovesOver(EliteSol_[0])) {
    return;
  }

  // Reject the new solution if it matches one of the elite solutions
  for (int idx=0; idx < EliteSol_.size(); ++idx) {
    if (x == EliteSol_[idx]) {
      return;
    }
  }

  // If we've made it to this point, we are going to add x to the elite set.

  // If the elite set is currently at capacity, remove the worst element
  // (maintaining the heap structure of the elite set). Also remove the
  // frequencies of the variables in the removed solution from EliteFreq_.
  if (EliteSol_.size() == R_) {
    const std::vector<int>& worst = EliteSol_[0].get_assignments();
    for (int j=0; j < N_; ++j) {
      if (worst[j] == 1) {
	--EliteFreq_[j];  // We're removing this solution, so decrement freq
      }
    }
    std::pop_heap(EliteSol_.begin(), EliteSol_.end(),
		  std::greater<Glover2010QUBOSolution>());
    EliteSol_.pop_back();
  }

  // Add x to the elite set, incrementing EliteFreq_ as appropriate
  const std::vector<int>& assignments = x.get_assignments();
  for (int j=0; j < N_; ++j) {
    if (assignments[j] == 1) {
      ++EliteFreq_[j];
    }
  }
  EliteSol_.push_back(x);
  std::push_heap(EliteSol_.begin(), EliteSol_.end(),
		 std::greater<Glover2010QUBOSolution>());
}

Glover2010::Glover2010(const QUBOInstance& qi, double runtime_limit,
		       bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  int R = 8;  // Maximum size of the memory EliteSol

  // Because Glover2010 is a population-based heuristic without an explicit
  // restart criterion, we will wrap it in random restarts. This is difficult
  // because it uses a runtime limit for its termination criterion. This limit
  // it set to require the same computational intensity as Palubeckis2006's
  // termination criterion, so we will use the exact same value as that
  // heuristic. Please see comments in palubeckis2006.cpp for details on how this
  // runtime limit was derived.
  double limit_raw;
  if (qi.get_size() < 2500) {
    limit_raw = 1.0 + 0.2690411 * qi.get_size();
  } else {
    limit_raw = 2490.8 - 1.635305 * qi.get_size() +
      0.0003633706 * qi.get_size() * qi.get_size();
  }
  double limit = limit_raw * 165.8905 / 1128.9522;

  while (true) {
    double start = Runtime();

    // Alg 1 Step 3: Initialize variables (EliteSol, r, and EliteFreq are all
    // maintained by Glover2010Elite)
    Glover2010Elite elite(qi, R);
    std::vector<int> FlipFreq;
    
    // Alg 1 Step 4: Initialize S0 to a random solution
    Glover2010QUBOSolution S0 = QUBOSolution::RandomSolution(qi, this);
    
    // Alg 1 Steps 5-14: Loop until R elite solutions are added (break was added
    // within the loop after the call to AddSolution().
    while (Runtime() - start <= limit) {
      // Alg 1 Step 6: Call tabu search on S0 (we'll keep the name S0 instead of
      // S*, which appears in the paper). FlipFreq will be modified by the tabu
      // search function as it operates.
      S0.TabuSearch(&FlipFreq);
      
      // Check termination criterion for algorithm
      if (!Report(S0)) {
        return;  // Out of time
      }
      
      // Alg 1 Steps 7-11: Insert S0 into the elite set unless it's already
      // there.
      elite.AddSolution(S0);
      
      // Break the loop if we just reached size R for the elite set (in the
      // pseudocode they do steps 12-13 before breaking, but this is repeated
      // again in steps 16-17, so it's just wasted effort.
      if (elite.size() == R) {
        break;
      }
      
      // Alg 1 Steps 12-13: Construct a new S0 as the perturbation operator
      // applied to a randomly selected elite solution
      S0 = Glover2010QUBOSolution::PerturbationOperator(elite.RandomSolution(),
                                                        elite.get_freq(),
                                                        elite.size(),
                                                        FlipFreq);
    }

    // Alg 1 Steps 15-24: Loop until runtime limit is reached.
    while (Runtime() - start <= limit) {
      // Alg 1 Steps 16-17: Construct a new S0 as the perturbation operator
      // applied to a randomly selected elite solution
      S0 = Glover2010QUBOSolution::PerturbationOperator(elite.RandomSolution(),
                                                        elite.get_freq(),
                                                        elite.size(),
                                                        FlipFreq);
      
      // Alg 1 Step 18: Call tabu search on S0 (we'll keep the name S0 instead of
      // S*, which appears in the paper). FlipFreq will be modified by the tabu
      // search function as it operates.
      S0.TabuSearch(&FlipFreq);
      
      // Check termination criterion for algorithm
      if (!Report(S0)) {
        return;  // Out of time
      }
      
      // Alg 1 Steps 19-23: Insert S0 into the elite set if it's not already there
      // and it improves on the worst solution in the elite set
      elite.AddSolution(S0);
    }
  }
}
