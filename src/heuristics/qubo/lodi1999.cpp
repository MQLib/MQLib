#include <math.h>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "heuristics/qubo/lodi1999.h"
#include "util/random.h"

Lodi1999MinRange::Lodi1999MinRange(const QUBOInstance& qi) :
  qi_(qi),
  CF_(qi.get_size(), -1),
  LB_(qi.get_lin()),
  UB_(qi.get_lin()) {
  // Compute LB_ and UB_ assuming no fixed indices
  for (auto iter = qi_.get_all_nonzero_begin();
         iter != qi_.get_all_nonzero_end(); ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double weight = iter->second;
    if (weight > 0) {
      UB_[i] += weight;
      UB_[j] += weight;
    } else {
      LB_[i] += weight;
      LB_[j] += weight;
    }
  }

  // Use the MinRangeInternal function to fix all the indices that we can,
  // updating LB_ and UB_. After this is run CF_ will be properly fixed.
  MinRangeInternal(&CF_, &LB_, &UB_, NULL);
}

void Lodi1999MinRange::MinRange(std::vector<int>* fixed,
                                ExtendedSolution* x) const {
  // Copy the LB_ and UB_ (the lower and upper bounds associated with the set
  // CF) into local vectors and pass those to MinRangeInternal.
  std::vector<double> local_LB(LB_);
  std::vector<double> local_UB(UB_);
  MinRangeInternal(fixed, &local_LB, &local_UB, x);
}

void Lodi1999MinRange::MinRangeInternal(std::vector<int>* fixed,
                                        std::vector<double>* LB,
                                        std::vector<double>* UB,
                                        ExtendedSolution* x) const {
  // Update bounds based on the passed fixed values
  for (int i=0; i < fixed->size(); ++i) {
    if ((*fixed)[i] != CF_[i]) {
      if (CF_[i] >= 0) {
        std::cout << "Mismatch with CF in MinRange" << std::endl;
        exit(1);
      } else {
        // Variable is fixed in this solution but not globally
        FixVariable(i, (*fixed)[i], LB, UB);
      }
    }
  }

  // Don't fix anything if we don't have a reference solution. Note that this
  // means we won't have any globally fixed indices. However, this is OK
  // because: 1) Those indices will be fixed immediately for all solutions, and
  // 2) When solving MAX-CUT (the current use of the code) there will never be
  // globally fixed variables.
  if (x == NULL) {
    return;
  }

  // Fix any additional variables until no more can be fixed
  bool updated = true;
  while (updated) {
    updated = false;
    for (int i=0; i < qi_.get_size(); ++i) {
      if ((*fixed)[i] >= 0) {
        continue;  // This variable is already fixed
      }
      if (x->ImprovingMove(-(*UB)[i])) {
        // Because we're maximizing, we will fix i to 0 instead of 1 as in paper
        (*fixed)[i] = 0;
        FixVariable(i, 0, LB, UB);
        updated = true;
      } else if (x->ImprovingMove((*LB)[i])) {
        // Because we're maximizing, we will fix i to 1 instead of 0 as in paper
        (*fixed)[i] = 1;
        FixVariable(i, 1, LB, UB);
        updated = true;
      }
    }
  }
}

void Lodi1999MinRange::FixVariable(int i, int val, std::vector<double>* LB,
                                   std::vector<double>* UB) const {
  // Loop through all values linked to i update their LB and UB
  for (auto iter = qi_.get_nonzero_begin(i); iter != qi_.get_nonzero_end(i);
       ++iter) {
    int j = iter->first;
    double weight = iter->second;
    if (weight > 0.0) {
      if (val == 0) {
        (*UB)[j] -= weight;
      } else {
        (*LB)[j] += weight;
      }
    } else {
      if (val == 0) {
        (*LB)[j] -= weight;
      } else {
        (*UB)[j] += weight;
      }
    }
  }
}

Lodi1999Solution::Lodi1999Solution(const QUBOSolution &x) :
QUBOSolution(x) {}

// phi iterations of the constructive phase (changing 0 to 1) and then
// destructive phase (changing 1 to 0) are performed. In each phase, we
// flip up to C values, at each inner iteration taking the one that will give us
// the biggest improvement in objective function.
void Lodi1999Solution::LS(const std::vector<int>& fixed) {
  // Parameters (see Lodi1999 constructor)
  const int phi = N_ < 200 ? 1 : 2;
  const int C = 50;
  for (int outer_iter=0; outer_iter < phi; ++outer_iter) {
    // Constructive phase
    for (int inner_iter=0; inner_iter < C; ++inner_iter) {
      // Find best diff_weight
      double max_dw = 0.0;
      int max_dw_i = -1;
      for (int i = 0; i < N_; i++) {
        if (fixed[i] >= 0 || assignments_[i] == 1) {
          continue;  // Ignore fixed values and 1s
        }
        if (diff_weights_[i] > max_dw) {
          max_dw = diff_weights_[i];
          max_dw_i = i;
        }
      }
      if (max_dw_i == -1 || !ImprovingMove(max_dw_i)) {
        break;  // No new indices to swap
      } else {
        UpdateCutValues(max_dw_i);
      }
    }

    // Destructive phase
    for (int inner_iter=0; inner_iter < C; ++inner_iter) {
      // Find best diff_weight
      double max_dw = 0.0;
      int max_dw_i = -1;
      for (int i = 0; i < N_; i++) {
        if (fixed[i] >= 0 || assignments_[i] == 0) {
          continue;  // Ignore fixed values and 0s
        }
        if (diff_weights_[i] > max_dw) {
          max_dw = diff_weights_[i];
          max_dw_i = i;
        }
      }
      if (max_dw_i == -1 || !ImprovingMove(max_dw_i)) {
        break;  // No new indices to swap
      } else {
        UpdateCutValues(max_dw_i);
      }
    }
  }
}

Lodi1999Solution::Lodi1999Solution(const Lodi1999Solution &father,
                                   const Lodi1999Solution &mother,
                                   const Lodi1999MinRange& mr) :
  QUBOSolution(father) {
  // PAPER: 4. generate a new solution by applying a cross-over operator to 
  //        the parents
  //        We temporarily fix the variables with the same value in P1 
  //        and P2 obtaining a child (SON) that has a reduced number of 
  //        unfixed variables and we apply to it the intensification actions
  std::vector<int> fixed(N_, -1);  // -1 if unfixed, 0/1 if both parents match
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] == mother.assignments_[i]) {
      fixed[i] = assignments_[i];
    }
  }

  // Use MinRange to update fixed to include all indices that should also
  // logically fixed to either 0 or 1.
  mr.MinRange(&fixed, this);

  // PAPER: 5. the new solution is modified by a mutation operator
  //        "Mutation operator: The mutation operator in genetic algorithms
  //         is intended to perform an action of diversification. In our
  //         algorithm we have two types of diversification:
  //         1. if at the end of the cross-over phase not all the
  //            variables have been fixed then the unfixed vari-
  //            ables of SON are randomly set;"
  for (int i = 0; i < N_; i++) {
    int value = (fixed[i] == -1) ? Random::RandInt(0, 1) : fixed[i];
    if (value != assignments_[i]) {
      UpdateCutValues(i);
    }
  }

  // PAPER: "Post-optimization: Between the phases of mu-
  //         tation and insertion, a post-optimization of the
  //         current solution SON is performed by applying
  //         algorithm LS"
  LS(fixed);
}

// Random solution followed by local search
Lodi1999Solution::Lodi1999Solution(const std::vector<int>& CF,
                                   const QUBOInstance& qi,
                                   QUBOHeuristic* heuristic) :
  QUBOSolution(qi, heuristic) {
  // Starting with empty solution, use index in CF if specified or otherwise
  // use random value.
  for (int i=0; i < N_; ++i) {
    if (CF[i] >= 0) {
      assignments_[i] = CF[i];
    } else {
      assignments_[i] = Random::RandInt(0, 1);
    }
  }
  PopulateFromAssignments();

  // Perform LS on the randomly generated solution
  LS(CF);
}

Lodi1999::Lodi1999(const QUBOInstance& qi, double runtime_limit,
		   bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {

  // READ ME VERY IMPORTANT
  // This paper is trying to solve a MINIMIZATION problem, so imp from
  // the paper (the vector of improvements from flipping each element) is
  // equal to the negative of diff_weights_ in our code.

  // NOTES: Class e and f are characterised by being larger than
  //        200 and 500 respectively, and have different parameters.
  //        Class a,b,c,d:
  //          C = 50, phi = 1(LS)
  //          SIZEmin = 50, SIZEmax = 80, R = 50, T = 3(EH)
  //        Class e:
  //          Increase phi to 2
  //        Class f:
  //          Increase phi to 2, SIZE_min = 60, R = 80
  const int SIZEmin = qi.get_size() < 500 ? 50 : 60;
  const int SIZEmax = 80;
  const int R = qi.get_size() < 500 ? 50 : 80;
  const int T = 3;

  std::vector<Lodi1999Solution> best_of_restarts;
  double global_best = -1e10;
  int no_global_best_improve_iterations = 0;

  // Preprocessing: compute all the globally fixed variables CF as well as the
  //                associated lower and upper bounds with this fixed set.
  Lodi1999MinRange mr(qi);
  int num_cf = 0;
  for (int i=0; i < qi.get_size(); ++i) {
    num_cf += mr.get_CF()[i] >= 0;
  }

  while (true) {
    // PAPER: 1: generate an initial population of different solutions, of which
    //           SIZE_min - g are randomly generated and improved through LS
    //           and the remaining g elements are the best solutions found in
    //           previous restarts performed by the algorithm.
    std::vector<Lodi1999Solution> population;
    for (int p = 0; p < SIZEmin - best_of_restarts.size(); p++) {
      population.push_back(Lodi1999Solution::RandomWithCF(mr.get_CF(), qi,
                                                          this));
    }
    for (int p = 0; p < best_of_restarts.size(); p++) {
      population.push_back(best_of_restarts[p]);
    }
    std::sort(population.begin(), population.end(),
              std::greater<Lodi1999Solution>());
    // PAPER: 2: repeat
    double best_weight = population[0].get_weight();
    int iterations_without_improvement = 0;
    while (true) {
      // PAPER: 3: select two solutions in the population, the parents
      //        "We consider the population sorted in non-decreasing order of
      //         the objective function value."
      //        "Note that by squaring y 1 and y 2 the probability of
      //         choosing the best solution in the population at
      //         each step becomes higher."
      // NOTES: We use decreasing instead of increasing order because we're
      //        maximizing instead of minimizing.
      // NOTES: It will always be sorted at this point
      double y1 = Random::RandDouble();
      double y2 = Random::RandDouble();
      int SIZE = population.size();
      int P1 = y1*y1 * SIZE + 1;
      int temp = y2*y2 * (SIZE - 1) + 1;
      int P2 = (temp < P1) ? temp : temp + 1;
      // PAPER: 4: generate a new solution by applying a cross-over operator
      //           to the parents
      //        5: the new solution is modified by a mutation operator
      assert(P1-1 >= 0);
      assert(P2-1 >= 0);
      assert(P1-1 < population.size());
      assert(P2-1 < population.size());
      Lodi1999Solution child = Lodi1999Solution::Child(population[P1-1],
                                                       population[P2-1], mr);
      // PAPER: 6: insert the new solution in the population and update the
      //           population
      // PAPER: "Insert and update: After having inserted SON in
      //         the current population, we sort the population and
      //         eliminate its worst solution if SIZE > SIZEmax ,
      //         where SIZEmax indicates the limit in the size of the
      //         population allowed."

      if (population.size() < SIZEmax ||
          child.ImprovesOver(population[population.size()-1])) {
        // We should insert the child into the sorted vector, removing the
        // worst solution if necessary.
        if (population.size() == SIZEmax) {
          population.pop_back();
        }
        std::vector<Lodi1999Solution>::iterator pos =
          std::lower_bound(population.begin(), population.end(),
                           child, std::greater<Lodi1999Solution>());
        population.insert(pos, child);
      }

      // PAPER: (ii) a prefixed time limit is reached
      assert(population.size() >= 1);
      if (!Report(population[0])) {
        break;
      }
      // PAPER: (i) a prefixed number of iterations is performed
      if (!population[0].ImprovesOver(best_weight)) {
        // Best didn't improve
        iterations_without_improvement++;
        if (iterations_without_improvement >= R) {
          // No improvement for too long
          break;
        }
      } else {
        best_weight = population[0].get_weight();
        iterations_without_improvement = 0;
      }
    } // end 2: repeat
    assert(population.size() >= 1);
    if (!Report(population[0])) {
      break;
    }
    // NOTES: We use the best of previous restarts in future populations
    best_of_restarts.push_back(population[0]);
    // PAPER: "after T restarts without improvement of the glob-
    //        ally best solution, we restart with a population
    //        containing no element of the previous restarts"
    if (!population[0].ImprovesOver(global_best)) {
      no_global_best_improve_iterations++;
      if (no_global_best_improve_iterations >= T) {
        no_global_best_improve_iterations = 0;
        best_of_restarts.clear();
      }
    } else {
      no_global_best_improve_iterations = 0;
      global_best = population[0].get_weight();
    }
  }
}
