#include <math.h>
#include <algorithm>
#include <iostream>
#include "heuristics/qubo/alkhamis1998.h"
#include "util/random.h"

void Alkhamis1998Solution::SA(double T_initial, int iteration) {
  // Parameters
  double SFACTOR = 0.3;
  double TFACTOR = 0.007;
  int ITER = std::max<int>(1, N_ * SFACTOR);
  double T_final = T_initial * TFACTOR;
  int ConsecFailureLimit = 10;
  double log_delta_plus_one = 0.09531018;  // log(1+delta) for delta = 0.1

  // Step 1: Setup state for the SA run
  int ConsecFailure = 0;
  double sd = DiffWeightStandardDeviation();  // We keep track of the sd
  double T = T_initial;  // Temperature

  // Step 2: Loop until frozen (either temp too low or too many consecutive
  //         failures).
  while (T > T_final && ConsecFailure < ConsecFailureLimit) {
    // Step 2.1: Test a random 1-swap ITER times
    double changed = false;  // Did we change in the ITER tries?
    Alkhamis1998Solution best(*this);
    for (int i=0; i < ITER; ++i) {
      // Step 2.1.1: Randomly select a variable to test
      int var = Random::RandInt(0, N_-1);

      // Steps 2.1.2-2.1.3: Test if move is accepted, and do it if it is
      if (NonDetrimentalMove(var) ||
          Random::RandDouble() < exp(diff_weights_[var] / T)) {
	UpdateCutValues(var);
	changed = true;

        // To avoid reporting in a tight loop, update the best solution if we have
        // obtained a new one (we'll report after this loop).
        if (ImprovesOver(best)) {
          best = *this;
        }
      }
    }

    // Report the best solution obtained during the loop
    if (!heuristic_->Report(best, iteration)) {
      return;
    }

    // Step 2.2.1: Update ConsecFailure and sd based on whether we changed this
    // iteration
    if (changed) {
      ConsecFailure = 0;
      sd = DiffWeightStandardDeviation();
    } else {
      ++ConsecFailure;
    }

    // Step 2.2.2: Update T according to formula (2)
    T /= (1.0 + T * log_delta_plus_one / 3.0 / sd);
  }
}

Alkhamis1998::Alkhamis1998(const QUBOInstance& qi, double runtime_limit,
			   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  double eta = 0.95;

  // Compute T_initial just once, since it's computationally intensive to compute
  // T_initial is computed by generating some undefined number "m" of trials.
  // We will set m=20 somewhat arbitrarily here.
  int m = 20;

  // The parameters to be computed (m1, m2, DeltaBarNegative) are all computed
  // with opposite sign from the paper, because the paper is minimizing but we
  // are maximizing.
  int m1 = 0;  // Variables with non-negative diff_weights_
  int m2 = 0;  // Variables with negative diff_weights_
  double DeltaBarNegativeSum = 0.0;  // Sum of negative diff_weights_ values

  for (int iter=0; iter < m; ++iter) {
    QUBOSolution x = QUBOSolution::RandomSolution(qi, this);
    const std::vector<double>& diff_weights = x.get_diff_weights();
    for (int i=0; i < diff_weights.size(); ++i) {
      if (x.NonDetrimentalMove(i)) {
	++m1;
      } else {
	++m2;
	DeltaBarNegativeSum += diff_weights[i];
      }
    }
  }

  // Finally, compute T_initial using formula (1)
  double T_initial =
    (-DeltaBarNegativeSum/m2) / log(m2 / (m2*eta - (1.0 - eta)*m1));

  // Random restart until termination criterion met
  for (int iter=0; QUBOHeuristic::Report(iter); ++iter) {
    // Steps 1.5-1.6: Randomly, generate a starting solution X and calculate f(X)
    Alkhamis1998Solution X(QUBOSolution::RandomSolution(qi, this));

    // Check termination criterion
    if (!Report(X, iter)) {
      return;
    }

    // Step 2: Run the simulated annealing procedure
    X.SA(T_initial, iter);
  }
}
