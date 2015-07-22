#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/qubo/lu2010.h"
#include "util/random.h"

Lu2010QUBOSolution::Lu2010QUBOSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

void Lu2010QUBOSolution::TabuSearch() {
  // Parameters (different values were given for "Random" and "SPP (Set
  // Partitioning Problem)); since random is QUBO we'll use those.
  int tt = (int)(N_ / 150);  // Tabu tenure constant
  int randTenure = 10;  // Upper bound on random component of tabu tenure
  int alpha = N_ * 5;  // Improvement cutoff for TS

  Lu2010QUBOSolution best = *this;  // Best solution yet encountered
  std::vector<int> TabuIter(N_, 0);  // Iteration when each var is non-tabu
  int numWithoutImprovement = 0;  // Moves w/o improving best solution
  for (int iter=0; numWithoutImprovement < alpha; ++iter) {
    std::vector<int> bests;  // All indices tied with best improvement
    double bestVal = -std::numeric_limits<double>::max();
    for (int i=0; i < N_; ++i) {
      // If var is either not on the tabu list or would improve on the best
      // solution to date, update "bests" and "bestVal"
      if (TabuIter[i] <= iter || ImprovesOverAfterMove(best, i)) {
        if (BaseSolution::ImprovesOver(weight_ + diff_weights_[i], bestVal)) {
	  bests.clear();
	  bests.push_back(i);
	  bestVal = weight_ + diff_weights_[i];
	} else if (!BaseSolution::ImprovesOver(bestVal,
                                               weight_ + diff_weights_[i])) {
	  bests.push_back(i);
	}
      }
    }

    // Randomly select index to flip from bests and flip it, updating the
    // tabu tenure for that index
    if (bests.size() > 0) {
      int idx = bests[Random::RandInt(0, bests.size() - 1)];
      UpdateCutValues(idx);
      TabuIter[idx] = iter + tt + Random::RandInt(1, randTenure) + 1;
    }

    // Check if we've improved on the best solution to date
    if (ImprovesOver(best)) {
      best = *this;
      numWithoutImprovement = 0;
    } else {
      ++numWithoutImprovement;
    }
  }

  // Copy over the best solution found during the tabu search to be the
  // current solution.
  *this = best;
}

Lu2010QUBOSolution::Lu2010QUBOSolution(const QUBOInstance& qi,
				       const Lu2010QUBOSolution& xi,
				       const Lu2010QUBOSolution& xj,
				       QUBOHeuristic *heuristic) :
  QUBOSolution(xi) {
  // Parameters
  double probDGPR = 0.5;  // Probability of using DG/PR instead of uniform
  if (Random::RandDouble() > probDGPR) {
    // Generate via uniform crossover (random value for mismatches)
    for (int i=0; i < N_; ++i) {
      if (xi.assignments_[i] != xj.assignments_[i] &&
	  Random::RandDouble() < 0.5) {
	// Mismatch and we decided to change the xi value
	UpdateCutValues(i);
      }
    }

  } else {
    // Generate via DG/PR.
    operator=(QUBOSolution(Lu2010PartialSolution::DGPR(qi, xi, xj, heuristic)));
  }
}

Lu2010PartialSolution::Lu2010PartialSolution(const QUBOInstance& qi,
					     const Lu2010QUBOSolution& xi,
					     const Lu2010QUBOSolution& xj,
					     QUBOHeuristic *heuristic) :
  QUBOPartialSolution(qi, heuristic) {
  // Identify the set NC (differing between xi and xj)
  std::vector<bool> isNC(N_, false);  // Does var differ between xi and xj?
  std::vector<int> NC;  // Indices that differ between xi and xj
  for (int i=0; i < N_; ++i) {
    if (xi.get_assignments()[i] != xj.get_assignments()[i]) {
      isNC[i] = true;
      NC.push_back(i);
      assignments_[i] = 0.5;
    } else {
      assignments_[i] = xi.get_assignments()[i];
    }
  }
  PopulateFromAssignments();

  // Iteratively assign the most promising move for each index.
  bool matchI = true;  // Are we matching to xi this round?
  while (NC.size() > 0) {
    // First, loop through the values in NC and pick the most promising
    int tomove = -1;  // The variable number to move
    int tomove_idx = -1;  // The index in NC of the variable to move
    double bestDiff = -std::numeric_limits<double>::max();
    for (int idx=0; idx < NC.size(); ++idx) {
      int i = NC[idx];
      double thisMove;
      if (matchI) {
	thisMove = xi.get_assignments()[i] ? diff1_[i] : diff0_[i];
      } else {
	thisMove = xj.get_assignments()[i] ? diff1_[i] : diff0_[i];
      }
      if (thisMove > bestDiff) {
	tomove = i;
	tomove_idx = idx;
	bestDiff = thisMove;
      }
    }
    
    // Assign tomove's variable value, updating isNC and NC.
    UpdateCutValues(tomove, matchI ? xi.get_assignments()[tomove] :
		    xj.get_assignments()[tomove]);
    isNC[tomove] = false;
    NC[tomove_idx] = NC[NC.size()-1];
    NC.resize(NC.size()-1);
    
    // Flip whether xi or xj is matched
    matchI = !matchI;
  }
}

Lu2010Population::Lu2010Population(int p, const QUBOInstance& qi,
				   QUBOHeuristic *heuristic) :
  p_(p),
  HD_(p*p, 0),
  NHD_(p*p, 0.0),
  min_NHD_(p, std::numeric_limits<double>::max()),
  VI_(qi.get_size(), 0.0) {
  // Parameters
  double phi = 0.2;  // Variable importance coefficient

  // Build a population by generating random solutions and then running tabu
  // search on each
  for (int ct=0; ct < p; ++ct) {
    P_.push_back(QUBOSolution::RandomSolution(qi, heuristic));
    P_[ct].TabuSearch();
    if (!heuristic->Report(P_[ct])) {
      return;  // Out of time
    }
  }

  // Compute the variable importance for each variable
  std::vector<double> interactionSums(qi.get_size(), 0.0);
  for (auto iter = qi.get_all_nonzero_begin(); iter != qi.get_all_nonzero_end();
       ++iter) {
    interactionSums[iter->first.first] += fabs(iter->second);
    interactionSums[iter->first.second] += fabs(iter->second);
  }
  for (int i=0; i < qi.get_size(); ++i) {
    VI_[i] = sqrt(fabs(qi.get_lin()[i]) + phi * interactionSums[i]);
  }

  // Compute the pairwise hamming distances (HD_), the pairwise non-hamming
  // distances (NHD_), min NHD to other elements (min_NHD_), and average hamming
  // distance (avg_HD_)
  double HD_sum = 0.0;
  for (int i=0; i < p; ++i) {
    for (int j=i+1; j < p; ++j) {
      std::vector<int> diffs;
      HD_[i*p + j] = P_[i].SymmetricDifference(P_[j], &diffs);
      HD_sum += HD_[i*p + j];
      for (int idx=0; idx < diffs.size(); ++idx) {
	NHD_[i*p + j] += VI_[diffs[idx]];
      }
      min_NHD_[i] = std::min<double>(min_NHD_[i], NHD_[i*p + j]);
      min_NHD_[j] = std::min<double>(min_NHD_[j], NHD_[i*p + j]);
    }
  }
  avg_HD_ = HD_sum * 2.0 / p / (p-1);
}

std::pair<const Lu2010QUBOSolution&, const Lu2010QUBOSolution&>
  Lu2010Population::RandomParents() const {

  // Check if all hamming distances are the same; if so an individual distance
  // can't exceed the average so we'll take any random pair of parents. Note
  // that this differs from the isConstant() function, which is checking if all
  // population elements are the same (aka they all have HD of 0).
  bool constHD = true;
  for (int i=0; i < p_; ++i) {
    for (int j=0; j < p_; ++j) {
      if (HD_[i*p_ + j] != HD_[0]) {
        constHD = false;
        break;
      }
    }
    if (!constHD) {
      break;
    }
  }

  while (true) {
    // Obtain random j != k
    int j = Random::RandInt(0, p_-1);
    int k = Random::RandInt(0, p_-1);
    while (j == k) {
      k = Random::RandInt(0, p_-1);
    }
    
    // Swap so j < k
    if (j > k) {
      int tmp = j;
      j = k;
      k = tmp;
    }

    // Return in success if the hamming distance of j and k exceeds the
    // average distance (or if all hamming distances are identical)
    if (HD_[j*p_ + k] > avg_HD_ || constHD) {
      return
	std::pair<const Lu2010QUBOSolution&, const Lu2010QUBOSolution&>(P_[j],
									P_[k]);
    }
  }
}

void Lu2010Population::UpdatePool(const Lu2010QUBOSolution& x0) {
  // Parameters
  double beta = 0.6;  // Goodness score coefficient
  double wp = 0.3;  // Probability of using new solution even if not good enough

  // First, we'll compute the non-Hamming distance (NHD) based on variable
  // importance from x0 to the other solutions in the population. We'll also
  // compute new_min_NHD, which is the min_NHD value for every current
  // element of the population, also including x0 in the computation.
  std::vector<int> x0_HD(p_, 0);
  std::vector<double> x0_NHD(p_, 0.0);
  std::vector<double> new_min_NHD = min_NHD_;
  // Min dist from x0 to an element in the population
  double x0_min_NHD = std::numeric_limits<double>::max();
  for (int i=0; i < p_; ++i) {
    std::vector<int> diffs;
    x0_HD[i] = x0.SymmetricDifference(P_[i], &diffs);
    for (int idx=0; idx < diffs.size(); ++idx) {
      x0_NHD[i] += VI_[diffs[idx]];
    }
    new_min_NHD[i] = std::min<double>(new_min_NHD[i], x0_NHD[i]);
    x0_min_NHD = std::min<double>(x0_min_NHD, x0_NHD[i]);
  }

  // Next, we'll determine the min and max values of the objective as well as
  // the distance to the population:
  double obj_min = x0.get_weight();
  double obj_max = x0.get_weight();
  double NHD_min = x0_min_NHD;
  double NHD_max = x0_min_NHD;
  for (int i=0; i < p_; ++i) {
    obj_min = std::min<double>(obj_min, P_[i].get_weight());
    obj_max = std::max<double>(obj_max, P_[i].get_weight());
    NHD_min = std::min<double>(NHD_min, new_min_NHD[i]);
    NHD_max = std::max<double>(NHD_max, new_min_NHD[i]);
  }

  // Now, compute the goodness score of x0 and each element in the population,
  // identifying w, the index of the worst-performing element in the population.
  double g_x0 = beta*((x0.get_weight() - obj_min) / (obj_max - obj_min + 1.0)) +
    (1-beta)*((x0_min_NHD - NHD_min) / (NHD_max - NHD_min + 1.0));
  int w = -1;
  double worst_score = std::numeric_limits<double>::max();
  for (int i=0; i < p_; ++i) {
    double g_xi = beta*((P_[i].get_weight()-obj_min)/(obj_max-obj_min + 1.0)) +
      (1-beta) * ((new_min_NHD[i] - NHD_min) / (NHD_max - NHD_min + 1.0));
    if (g_xi < worst_score) {
      worst_score = g_xi;
      w = i;
    }
  }

  // If x0 is selected to be added to the population, update class values.
  if (!BaseSolution::ImprovesOver(worst_score, g_x0) ||
      Random::RandDouble() < wp) {
    // Insert x0 into position w, updating P_, HD_, and NHD_.
    P_[w] = x0;
    for (int i=0; i < w; ++i) {
      HD_[i*p_ + w] = x0_HD[i];
      NHD_[i*p_ + w] = x0_NHD[i];
    }
    for (int i=w+1; i < p_; ++i) {
      HD_[w*p_ + i] = x0_HD[i];
      NHD_[w*p_ + i] = x0_NHD[i];
    }

    // Compute avg_HD_ and min_NHD_ for the new population
    for (int i=0; i < p_; ++i) {
      min_NHD_[i] = std::numeric_limits<double>::max();
    }
    double HD_sum = 0.0;
    for (int i=0; i < p_; ++i) {
      for (int j=i+1; j < p_; ++j) {
	HD_sum += HD_[i*p_ + j];
	min_NHD_[i] = std::min<double>(min_NHD_[i], NHD_[i*p_ + j]);
	min_NHD_[j] = std::min<double>(min_NHD_[i], NHD_[i*p_ + j]);
      }
    }
    avg_HD_ = HD_sum * 2.0 / p_ / (p_-1);
  }
}

Lu2010::Lu2010(const QUBOInstance& qi, double runtime_limit, bool validation,
               QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  double p = 20;  // Population size
  
  // Outer iteration: once the population converges to being constant there's
  // little point in continuing, so we'll instead restart whenever that occurs.
  while (true) {
    
    // Build the initial population, by randomly generating p solutions and
    // performing tabu search on each
    Lu2010Population P(p, qi, this);
    // If we hit the termination criterion while generating the population, just
    // exit.
    if (!QUBOHeuristic::Report()) {
      break;
    }
    
    // Alg 1, Step 8 and 16: Repeat until termination criterion met. There's no
    // point in continuing after all solutions in P are the same, so we'll use
    // that as the inner termination criterion.
    while (!P.isConstant()) {
      // Alg 1, Step 9: Randomly select individuals j and k for combination,
      // limiting to pairs that are above-average diversity from one another in P.
      std::pair<const Lu2010QUBOSolution&, const Lu2010QUBOSolution&> parents =
        P.RandomParents();

      // Alg 1, Step 10: Build x0 from xj and xk
      Lu2010QUBOSolution x0 = Lu2010QUBOSolution::Combine(qi, parents.first,
                                                          parents.second, this);

      // Alg 1, Step 11: Perform tabu search on x0
      x0.TabuSearch();

      // Alg 1, Steps 12-14: Check if new best encountered (we'll also check our
      // termination criterion here)
      if (!Report(x0)) {
        return;  // Time is up
      }

      // Alg 1, Step 15: Run Pool_Updating procedure with new solution
      P.UpdatePool(x0);
    }
  }
}
