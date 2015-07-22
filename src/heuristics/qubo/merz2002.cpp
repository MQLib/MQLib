#include <algorithm>
#include <iostream>
#include "heuristics/qubo/merz2002.h"
#include "util/random.h"

Merz2002PartialSolution::Merz2002PartialSolution(const QUBOInstance& qi,
						 QUBOHeuristic *heuristic) :
  QUBOPartialSolution(qi, heuristic) {
  // Initialize current solution to all assignments being 0.5.
  assignments_.assign(N_, 0.5);
  PopulateFromAssignments();

  // Select random variable and randomly select value to either 0 or 1.
  UpdateCutValues(Random::RandInt(0, N_-1), Random::RandInt(0, 1));

  // Iterate until no variables have fractional values
  while (num_frac_ > 0) {
    // Identify the fractional vars with the best gain to 0 or 1
    double best0val = -std::numeric_limits<double>::max();
    int best0pos = -1;
    double best1val = -std::numeric_limits<double>::max();
    int best1pos = -1;
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] == 0 || assignments_[i] == 1.0) {
	continue;  // Not a fractional solution
      }
      if (diff0_[i] > best0val) {
	best0val = diff0_[i];
	best0pos = i;
      }
      if (diff1_[i] > best1val) {
	best1val = diff1_[i];
	best1pos = i;
      }
    }

    // Randomly select which to update, based on their weights. ** Because
    // merz2002 does not comment on what to do if one or both of these values
    // are negative, we'll do something natural (if both are negative, pick
    // randomly; if one is negative, update the positive one). **
    double flip1;  // We will flip the best choice for 1
    if ((best0val > 0.0 && best1val > 0.0) ||
        (best0val == 0.0 && best1val == 0.0)) {
      flip1 = (Random::RandInt(0, 1) == 1);
    } else if (best0val > 0.0) {
      flip1 = false;
    } else if (best1val > 0.0) {
      flip1 = true;
    } else {
      flip1 = (Random::RandDouble() > (best1val / (best0val + best1val)));
    }

    // Assign the selected variable to its non-fractional value
    if (flip1) {
      UpdateCutValues(best1pos, 1);
    } else {
      UpdateCutValues(best0pos, 0);
    }
  }
}

void Merz2002QUBOSolution::KOpt() {
  // best is the best solution found in the k-opt. At the end of each iteration,
  // we will copy over this best solution to be the current solution.
  Merz2002QUBOSolution best(*this);

  // Loop until no improvement during loop (we will break there)
  bool improved = true;
  while (improved) {
    // Initialize iteration variables (in terminology of Figure 6, improved
    // becomes true when Gmax > 0.
    improved = false;
    std::vector<bool> inC(N_, true);  // Is a variable in set C (not flipped)?
    for (int iter=0; iter < N_; ++iter) {
      // Identify the best 1-flip to perform
      int bestPos = -1;
      double bestDW = -std::numeric_limits<double>::max();
      for (int i=0; i < N_; ++i) {
	if (inC[i] && diff_weights_[i] > bestDW) {
	  bestPos = i;
	  bestDW = diff_weights_[i];
	}
      }
      
      // Perform the 1-flip and update inC
      UpdateCutValues(bestPos);
      inC[bestPos] = false;
      
      // If we have improved over best with this flip, update best and improved.
      if (ImprovesOver(best)) {
	best = *this;
	improved = true;
      }
    }

    // Copy over the best solution to be the current solution
    operator=(best);
  }
}

Merz2002::Merz2002(const QUBOInstance& qi, double runtime_limit, bool validation,
		   QUBOCallback *qc, bool greedy, LStype type) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Add random restarts around procedure until out of time
  while (1) {
    // Step 1: Construct procedure (either random or greedy)
    Merz2002QUBOSolution x =
      greedy ? Merz2002PartialSolution::RandomizedGreedy(qi, this) :
      QUBOSolution::RandomSolution(qi, this);

    // Step 2: Local search (either none, best 1-opt, or k-opt)
    switch (type) {
    case ONEOPT:
      x.AllBest1Swap();
      break;
    case KOPT:
      x.KOpt();
      break;
    default:
      break;
    }

    // Check termination criterion with new solution
    if (!Report(x)) {
      break;
    }
  }
}
