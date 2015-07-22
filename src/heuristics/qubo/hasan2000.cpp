#include <algorithm>
#include <iostream>
#include <limits>
#include "heuristics/qubo/hasan2000.h"
#include "util/random.h"

Hasan2000Solution::Hasan2000Solution(const Hasan2000Solution& x1,
				     const Hasan2000Solution& x2) :
  QUBOSolution(x1) {
  // For a randomly selected crossover point k, we will take positions k+1
  // and k+2 from x2 and the rest from x1 (assuming crossingLength 2)
  int k = Random::RandInt(0, N_-1);
  int limit = std::min<int>(N_-1, k+2);
  for (int i=k+1; i <= limit; ++i) {
    if (x2.assignments_[i] != x1.assignments_[i]) {
      UpdateCutValues(i);
    }
  }
}

void Hasan2000Solution::Mutate() {
  // Parameters
  double mutationRate = 1.0 / N_;

  // Mutate each node with probability mutationRate
  for (int i=0; i < N_; ++i) {
    if (Random::RandDouble() < mutationRate) {
      UpdateCutValues(i);
    }
  }
}

Hasan2000Elite::Hasan2000Elite(const QUBOInstance& qi, int POP,
			       QUBOHeuristic *heuristic) :
  POP_(POP) {
  // Add POP random solutions, and perform local search on each
  for (int iter=0; iter < POP; ++iter) {
    P_.push_back(QUBOSolution::RandomSolution(qi, heuristic));
    P_[iter].AllBest1Swap();
  }
}

bool Hasan2000Elite::Update(const Hasan2000Solution& x) {
  // Exit if the new solution matches any existing solution in the population.
  // Also use this loop to identify all elite solutions that are worse than
  // the passed solution.
  std::vector<int> worse;
  for (int i=0; i < POP_; ++i) {
    if (x == P_[i]) {
      return false;  // Duplicated
    }
    if (x.ImprovesOver(P_[i])) {
      worse.push_back(i);
    }
  }

  // The paper says "the child replaces a poorer member of the population". We
  // take this to mean the child replaces a randomly selected element worse than
  // it.
  if (worse.size() > 0) {
    int position = worse[Random::RandInt(0, worse.size()-1)];
    P_[position] = x;
  }
  return true;  // Not duplicated
}

Hasan2000GA::Hasan2000GA(const QUBOInstance& qi, double runtime_limit,
			 bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Parameters
  int iterBeforeMutation = std::max<int>(1, 0.2 * qi.get_size());
  int numTournament = 4;
  int POP = 100;  // population size
  int nonDuplicateLimit = 20000;

  // Wrap the whole procedure in a random restart, restarting every time we have
  // seen at least nonDuplicateLimit number of non-duplicated children solutions.
  while (true) {
    // Generate the initial population (Section 4.1.1)
    Hasan2000Elite P(qi, POP, this);
    
    // Continue until we've seen a specified number of non-duplicate children
    // generated.
    int numCrossover = 0;
    int numNonDuplicate = 0;
    while (numNonDuplicate < nonDuplicateLimit) {
      // Section 4.1.2: Run numTournament tournaments, getting a child from each.
      // Binary tournaments are used, which is just rand selection of 2 parents.
      std::vector<Hasan2000Solution> children;
      for (int iter = 0; iter < numTournament; ++iter) {
        double i = Random::RandInt(0, POP-1);
        double j;
        while ((j = Random::RandInt(0, POP-1)) == i);
        
        // Section 4.1.3: One-point crossover with a specified crossing length.
        // If this child is selected for mutation (according to number of
        // crossovers so far), then do the mutation -- see Section 4.1.4.
        children.push_back(Hasan2000Solution::Crossover(P.get_P()[i],
                                                        P.get_P()[j]));
        ++numCrossover;
        if (numCrossover % iterBeforeMutation == 0) {
          children[children.size() - 1].Mutate();
        }
        children.push_back(Hasan2000Solution::Crossover(P.get_P()[j],
                                                        P.get_P()[i]));
        ++numCrossover;
        if (numCrossover % iterBeforeMutation == 0) {
          children[children.size() - 1].Mutate();
        }
      }
      
      // The paper doesn't specify it, but almost certainly they want to do a
      // local search here -- it's what they do to the initial solutions and it
      // makes the algorithm competitive.
      for (int i=0; i < children.size(); ++i) {
        children[i].AllBest1Swap();
      }
      
      // Identify the best child, adding each to the population if it doesn't
      // duplicate a solution currently in the population.
      int bestIdx = -1;
      double bestWeight = -std::numeric_limits<double>::max();
      for (int i=0; i < children.size(); ++i) {
        if (children[i].get_weight() > bestWeight) {
          bestIdx = i;
          bestWeight = children[i].get_weight();
        }
        if (P.Update(children[i])) {
          ++numNonDuplicate;
        }
      }
      
      // Exit if out of time
      if (!Report(children[bestIdx])) {
        return;
      }
    }
  }
}

void Hasan2000Solution::TS() {
  // Parameter: length of time spent of tabu list after move; the paper used
  // N_, but that doesn't make any sense (this means after a while only a
  // single node is considered at each iteration). Therefore, we use N_/2, which
  // it the value proposed for the BA strategy. Of course, this has much
  // better results.
  int TLs = N_/2;
  int MAXI = 4 * N_;

  // non-tabu if the iteration count exceeds the TABL value
  std::vector<int> TABL(N_, -1);
  // Keep track of best solution so we don't need to call Report in a tight loop
  Hasan2000Solution best(*this);
  int nonImproving = 0;  // Number of non-improving iterations
  for (int t=0; nonImproving < MAXI; ++t) {
    // Determine the best move to perform, limiting to non-tabu moves and
    // moves that meet the aspiration criterion (improving the best ever sln).
    // We use the FA strategy, so if a move improves the current objective we
    // select it without checking further.
    int selectedMove = -1;
    double selectedDW = -std::numeric_limits<double>::max();
    for (int i=0; i < N_; ++i) {
      if (t <= TABL[i]) {
	// i is tabu
	if (ImprovesOverAfterMove(best, i)) {
	  selectedMove = i;
	  break;
	}
      } else {
	// i is not tabu
        if (ImprovingMove(i)) {
	  selectedMove = i;
	  break;
	} else if (diff_weights_[i] > selectedDW) {
	  selectedMove = i;
	  selectedDW = diff_weights_[i];
	}	
      }
    }

    // Perform our selected move, updating the tabu information for this var
    if (selectedMove >= 0) {
      UpdateCutValues(selectedMove);
      TABL[selectedMove] = t + TLs;
    }

    // Update the best solution but only report periodically to avoid reporting
    // in a tight loop.
    if (ImprovesOver(best)) {
      best = *this;
      nonImproving = 0;
    } else {
      ++nonImproving;
    }
    if (t % 1000 == 0 && !heuristic_->Report(best)) {
      return;
    }
  }
}

Hasan2000TS::Hasan2000TS(const QUBOInstance& qi, double runtime_limit,
			 bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Create a random solution and run the tabu search on it. Random restart added
  // to prevent premature convergence.
  while (true) {
    Hasan2000Solution x = QUBOSolution::RandomSolution(qi, this);
    x.TS();

    if (!Report(x)) {
      break;
    }
  }
}
