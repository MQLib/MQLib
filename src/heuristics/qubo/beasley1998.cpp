#include <math.h>
#include <algorithm>
#include <iostream>
#include "heuristics/qubo/beasley1998.h"
#include "util/random.h"

Beasley1998Solution::Beasley1998Solution(const QUBOSolution &x) :
QUBOSolution(x) {}

void Beasley1998Solution::SA(double T) {
  // PAPER: randomly select a variable k
  int k = Random::RandInt(0, N_ - 1);
  // PAPER: check if new solution better than current solution
  // NOTES: This captures two cases that were separated in the papers (where it
  //        improves over the best solution ever and when it doesn't do that but
  //        does improve over the current solution). We combine those steps here
  //        and handle updating the best solution ever in the main loop.
  if (ImprovingMove(k)) {
    // NOTES: If it's better, just take it
    UpdateCutValues(k);
  } else {
    // PAPER: V** is the solution value associated with the current solution
    // NOTES: By current solution, its means the incumbent before any swap.
    //        V is the value of a solution after the swap has been made for k
    //        We never explicitly form this.
    //        Therefore:
    //          V** = weight_
    //          V   = weight_ + diff_weights_[k]
    //        So:
    //          V** - V = diff_weights_[k]
    if (Random::RandDouble() < exp(diff_weights_[k] / T)) {
      // NOTES: Accept worse solution
      UpdateCutValues(k);
    }
  }
}

void Beasley1998Solution::LocalSearch(int &t) {
  // NOTES: This is essentially AllFirst1Swap, with the catch that
  //        we are increasing an iteration counter that is shared 
  //        with the main algorithm in Tabu Search.
  bool improved = true;
  while (improved) {
    improved = false;
    for (int k = 0; k < N_; k++) {
      t++;
      if (ImprovingMove(k)) {
        UpdateCutValues(k);
        improved = true;
      }
    }
  }
}

// L is the first iter where it's OK to use each vertex, iter is the current
// iteration, vStar is the best solution found in this run of the algorithm,
// and t is a reference to the operation count, which is used to terminate the
// procedure.
int  Beasley1998Solution::TS(std::vector<int> &L, int iter, double vStar,
			     int &t) {
  int K = -1;
  
  // PAPER: Initialise best neighbour value
  double V_starstar = -1e10;
  // PAPER: examine all non-tabu variables
  for (int k = 0; k < N_; k++) {
    if (L[k] > iter) continue;
    t++;

    // NOTES: V is what the objective value would be if we swapped k. We never
    //        explictly form this solution.
    double V = weight_ + diff_weights_[k];

    // NOTES: If flipping k make this the best solution we have encountered,
    //        then take the move and local search
    if (ImprovesOver(V, vStar)) {
      // NOTES: now we actually form the swapped solution
      UpdateCutValues(k);
      // PAPER: apply the local search procedure 
      LocalSearch(t);
      // PAPER: go to done (tabu update handled in main loop)
      return k;
    }

    // NOTES: get here if doesn't improve on the best solution ever if we flip k. 
    // PAPER: check for improved neighbouring solution
    if (V > V_starstar) {
      K = k;
      V_starstar = V;
    }
  }

  if (K != -1) {
    // PAPER: make the move for chosen variable K
    // NOTES: and there is a move to make ()
    UpdateCutValues(K);
  }
  
  // PAPER: done:
  // NOTES: do tabu update in main loop
  return K;
}

Beasley1998SA::Beasley1998SA(const QUBOInstance& qi, double runtime_limit,
			     bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {

  // NOTES: Random restart until we hit the termination criterion
  for (int iter=0; QUBOHeuristic::Report(iter); ++iter) {
    // PAPER: randomise the starting solution
    Beasley1998Solution sol(QUBOSolution::RandomSolution(qi, this));
    // NOTES: We also need to store the best solution we found so far
    Beasley1998Solution best_sol(sol);

    // PAPER: initialise SA parameters
    double T = qi.get_size();
    const double alpha = 0.995;
    // PAPER: set maximum number of iterations
    const int T_star = std::max(500000, 5000 * qi.get_size());
    // PAPER: initialise iteration counter
    int t = 0;
    // PAPER: T* iterations in all
    while (t < T_star) {
      // PAPER: increment iteration counter
      t++;
      // NOTES: All the SA step is done as function of solution
      sol.SA(T);
      // PAPER: reduce temperature
      T = alpha * T;
      // NOTES: Need to update best ever, which we'll apply local search on
      //        at the end of the procedure
      if (sol.ImprovesOver(best_sol)) {
        best_sol = sol;
      }

      // To avoid calling Report in a tight loop, only check every 10000
      // iterations
      if (t % 10000 == 0) {
        if (!Report(best_sol, iter)) {
          return;  // Exit if termination criterion met
        }
      }
    }

    // PAPER: apply the local search procedure
    best_sol.LocalSearch(t);
    if (!Report(best_sol, iter))
      return;  // Exit if termination criterion met
  }
}

Beasley1998TS::Beasley1998TS(const QUBOInstance& qi, double runtime_limit,
			     bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {

  // NOTES: Restart until out of time unless it's a validation run, in which case
  // we just perform one run until we hit the iteration limit T_star.
  for (int iter=0; QUBOHeuristic::Report(iter); ++iter) {
    // PAPER: set the starting solution X_i = 0  i =1...n
    //Beasley1998Solution sol(QUBOSolution::UniformSolution(qi, 0, this));
    // NOTES: This was changed to random starts so the algorithm can 
    //        run forever
    Beasley1998Solution sol(QUBOSolution::RandomSolution(qi, this));
    // PAPER: initialise best solution value
    double vStar = sol.get_weight();
    // PAPER: initialise tabu values (for efficiency, store the first iteration
    //        for which this vertex will be non-tabu)
    std::vector<int> L(qi.get_size(), 0);
    // PAPER: set tabu tenure
    const int L_star = std::min(20, qi.get_size()/4);
    // PAPER: set maximum number of iterations
    const int T_star = std::max(500000, 5000 * qi.get_size());
    // PAPER: intialise iteration counter
    int t = 0;
    // PAPER: T* iterations in all (count inner iterations for efficient tabu
    //        search implementation)
    for (int inner_iter = 0; t < T_star; ++inner_iter) {
      int K = sol.TS(L, inner_iter, vStar, t);

      // PAPER: record improved solution (moved out of local search step)
      if (sol.get_weight() > vStar) {
	vStar = sol.get_weight();
      }

      // PAPER: tabu the chosen variable (L stores next iter where it's OK to use)
      if (K != -1)
        L[K] = inner_iter + L_star + 1;

      if (!Report(sol, iter))
        return;  // Exit if termination criterion met
    }
    if (!Report(sol, iter))
      return;  // Exit if termination criterion met
  }
}
