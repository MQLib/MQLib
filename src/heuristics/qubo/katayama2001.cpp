#include <math.h>
#include <algorithm>
#include <iostream>
#include "heuristics/qubo/katayama2001.h"
#include "util/random.h"

Katayama2001Solution::Katayama2001Solution(const QUBOSolution &x) :
QUBOSolution(x) {}

bool Katayama2001Solution::SASwap(int k, double T) {
  if (ImprovingMove(k)) {
    // PAPER:   ...then set Counter = 0...    
    // PAPER:   ...x_k = 1 - x_k (and update all gain g_i)...
    UpdateCutValues(k);
    // PAPER:   ...ignore rest
    return true;
  } else {
    // PAPER: SA: 3.3.3 Otherwise, set x_k = 1 - x_k with
    //                  probability exp(g_k/T)
    if (Random::RandDouble() < exp(diff_weights_[k] / T)) {
      UpdateCutValues(k);
    }
    return false;
  }
}

Katayama2001::Katayama2001(const QUBOInstance& qi, double runtime_limit,
			   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {

  std::vector<int> RP(qi.get_size(), 0);
  for (int i = 0; i < qi.get_size(); i++)
    RP[i] = i;

  // PAPER: Fig. 4. SA algorithm with reannealing process.
  // PAPER: 1 Select a ...
  double T_INIT = 0.3 * static_cast<double>(qi.get_size());
  const double T_FACTOR = 0.99;
  const int TERM_COUNT = 10;
  const int SA_COUNT = 2;
  const double START_T_FACTOR = 0.8;

  // PAPER: 2 Generate an initial random solution x_best
  // NOTES: Loop until termination criterion met (only one loop for valdiation)
  for (int iter=0; QUBOHeuristic::Report(iter); ++iter) {
    Katayama2001Solution x_best = QUBOSolution::RandomSolution(qi, this);
    if (!Report(x_best, iter))
      return;
    
    // PAPER: 3 Perform the following SACount times
    for (int outer_iter = 0; outer_iter < SA_COUNT; outer_iter++) {
      // PAPER: 3.1 x = x_best
      // PAPER: 3.2 x_best = SA(x,T_init,TFactor,TermCount)
      // PAPER: SA: 1 Set T = T_init, x_best = x, Counter = 0
      double T = T_INIT;
      int counter = 0;
      // PAPER: SA: 2 Calculate gains g[i] of x for all i in {1,...,n}
      // NOTES: This is already done by the code inside QUBOInstance
      // PAPER: SA: 3 Do the following until TermCount < Counter
      while (counter < TERM_COUNT) {
        // PAPER: SA: 3.1 Counter = Counter + 1
        counter++;
        // PAPER: SA: 3.2 Generate a random permutation RP[] ranging
        //                from 1 to n
        std::random_shuffle(RP.begin(), RP.end());
        // PAPER: SA: 3.3 For j = 1 to n
        for (int j = 0; j < qi.get_size(); j++) {
          // PAPER: SA: 3.3.1 k = RP[j];
          int k = RP[j];
          // PAPER: SA: 3.3.2 
          // PAPER:   If g_k > 0...
          if (x_best.SASwap(k, T)) {
            counter = 0;
          }
          // PAPER: SA: 3.4 T = TFactor x T
          T *= T_FACTOR;
        }
        // Check we don't run out of time inside SA
        if (!Report(x_best, iter))
          return;
      }  // End loop of Func SA(...)
      // PAPER: SA: 4 Return x_best
      // NOTES: End of procedure Sim-Anneal line 3.2 / func SA(...)
      // PAPER: 3.3 Set T_init = StartTFactor x T_init
      T_INIT *= START_T_FACTOR;
      // PAPER: 4 Return x_best
      if (!Report(x_best, iter))
        return;  // Exit if out of time and not doing validation run
    }
  }
}
