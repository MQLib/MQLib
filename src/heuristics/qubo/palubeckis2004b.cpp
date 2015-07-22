#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include "heuristics/qubo/palubeckis2004b.h"
#include "util/random.h"

Palubeckis2004bSolution::Palubeckis2004bSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

void Palubeckis2004bSolution::STSMod(double *best_objective, int zmax,
				     std::vector<Palubeckis2004bSolution>* B,
				     int mStar, double Delta,
				     Palubeckis2004bSolution* best,
				     bool reportBest) {
  // Step 0: Parameters for STS method
  int T = std::min(20, N_/4);  // Tabu period for flipped nodes

  // Step 1: Initialize tabu list and operation counter
  // Tabu list -- for efficiency store next iteration where allowed
  std::vector<int> TabuList(N_, 0);
  int z = 0;  // Operation counter
  for (int iter=0; z < zmax; ++iter) {
    // Step 2: Initialize some iteration variables
    double V = -std::numeric_limits<double>::max();  // Best possible change
    int r = -1;  // Best possible index
    int rho = 0;  // Have we found a best ever solution?
    for (int k=0; k < N_; ++k) {
      // Step 3.1: Don't process things on the tabu list
      if (TabuList[k] > iter) {
	continue;
      }

      // Step 3.2: Increment operation counter
      ++z;

      // Step 3.3: Check if this move would make a best ever solution
      if (ImprovesOverAfterMove(*best_objective, k)) {
	r = k;
	rho = 1;
	break;
      }

      // Step 3.4: Check if this move is the best found this iteration
      if (diff_weights_[k] > V) {
	V = diff_weights_[k];
	r = k;
      }
    }

    // Step 4: Take the selected move
    UpdateCutValues(r);

    // Step 5: If we just found the best ever solution, perform local search and
    // update the best ever record.
    if (rho) {
      LocalSearch(&z);
      *best_objective = weight_;  // New best objective ever
      if (best != NULL) {
	*best = *this;
      }
      if (reportBest && !heuristic_->Report(*this)) {
	return;  // Out of time
      }
    }

    // Adaptive Memory update step
    if (B != NULL) {
      if (B->size() == 0) {
	B->push_back(QUBOSolution(*this));
      } else {
	// Loop through B, checking if our current solution matches a solution
	// already there. Also extract the worst_pos and worst_weight, the
	// index and weight of the worst element of B.
	bool matches = false;
	int worst_pos = -1;
	double worst_weight = std::numeric_limits<double>::max();
	for (int ct=0; ct < B->size(); ++ct) {
	  if ((*B)[ct] == *this) {
	    matches = true;
	    break;
	  }
	  if ((*B)[ct].weight_ < worst_weight) {
	    worst_pos = ct;
	    worst_weight = (*B)[ct].weight_;
	  }
	}

	// Determine if we want to include this element
	bool include;
	if (!matches) {
	  if (ImprovesOver(worst_weight)) {
	    include = true;  // Always include if improving on some elt of B
	  } else if (BaseSolution::ImprovesOver(worst_weight, weight_) &&
                     B->size() >= mStar) {
	    include = false;  // Exclude if B is full and solution worse than all
	  } else {
	    // Check diversity condition
	    include = true;
	    for (auto iter=B->begin(); iter != B->end(); ++iter) {
	      if (SymmetricDifference(*iter) <= Delta * N_) {
		include = false;
		break;
	      }
	    }
	  }
	} else {
	  include = false;
	}

	// Actually add the current solution to the B set
	if (include) {
	  if (B->size() < mStar) {
	    B->push_back(QUBOSolution(*this));
	  } else {
	    (*B)[worst_pos] = QUBOSolution(*this);
	  }
	}
      }
    }

    // Step 6: Update the tabu list
    TabuList[r] = iter + T + 1;
  }
}

void Palubeckis2004bSolution::LocalSearch(int *z) {
  while (1) {
    // Step 1: Initialize rho (whether we've improved this iteration
    int rho = 0;

    for (int r=0; r < N_; ++r) {
      // Step 2.1: Increment operation counter
      ++(*z);

      // Step 2.2: Perform this 1-move if it's improving
      if (ImprovingMove(r)) {
	UpdateCutValues(r);
	rho = 1;
      }
    }

    // Step 3: Repeat only if we made an improvement this iteration
    if (!rho) {
      break;
    }
  }
}

Palubeckis2004bMST1::Palubeckis2004bMST1(const QUBOInstance& qi,
					 double runtime_limit, bool validation,
                                         QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z;
  if (qi.get_size() >= 2500) {
    Z = 15000;
  } else if (qi.get_size() >= 1000) {
    Z = 13000;
  } else {
    Z = 3500;
  }
  int zmax = std::max(500000, qi.get_size() * Z);
  
  double best_objective = 0.0;
  bool first = true;
  while (true) {
    // Step 2: Randomly assign variables as 1 or 0. If this is the first run,
    // assign it as the best objective.
    Palubeckis2004bSolution solution = QUBOSolution::RandomSolution(qi, this);
    if (first) {
      first = false;
      best_objective = solution.get_weight();
    }
    
    // Step 3: Run STS
    solution.STS(&best_objective, zmax);
    
    // Check if out of time
    if (!Report(solution)) {
      break;
    }
  }
}

void Palubeckis2004bSolution::SelectVariables(int n_prime,
					      std::vector<int>* I_star) const {
  // Parameters
  double lambda = 5000.0;

  // Can't select more than all nodes
  n_prime = std::min<int>(n_prime, N_);
  // First value is d_i, other is whether we have already added i to I_star
  std::vector<std::pair<double, bool> > d;

  // Step 1: I* is empty, d values set to diff_weights_ values.
  I_star->clear();
  for (int ct=0; ct < N_; ++ct) {
    d.push_back(std::pair<double, bool>(diff_weights_[ct], false));
  }

  // Loop until n_prime values added to I_star
  for (int ct=0; ct < n_prime; ++ct) {
    // Step 2: compute e values from d values. First, grab d_min and d_max:
    double d_min = std::numeric_limits<double>::max();
    double d_max = -std::numeric_limits<double>::max();
    for (int ct=0; ct < N_; ++ct) {
      if (!d[ct].second) {
	d_max = std::max<double>(d_max, d[ct].first);
	d_min = std::min<double>(d_min, d[ct].first);
      }
    }

    // Now, actually compute the e values.
    std::vector<double> e_value;
    std::vector<int> e_index;
    if (d_max - d_min < 1e-6) {
      for (int i=0; i < N_; ++i) {
	if (!d[i].second) {
	  e_value.push_back(1.0);
	  e_index.push_back(i);
	}
      }
    } else {
      for (int i=0; i < N_; ++i) {
	if (!d[i].second) {
	  if (d[i].first <= 0.0 && d_min < 0.0) {
	    e_value.push_back(1.0 - d[i].first / d_min);
	  } else if (d[i].first == 0.0 && d_min == 0.0) {
	    e_value.push_back(0.0);
	  } else {
            // d[i].first positive (and hence d_max positive)
	    e_value.push_back(1.0 + lambda * d[i].first / d_max);
	  }
	  e_index.push_back(i);
	}
      }
    }

    // Step 3: select k by roulette wheel selection
    int k = e_index[Random::RouletteWheel(e_value)];

    // Step 4: Add k to I_star. Update d based on selected k value. Note that
    // c'_{ik} = c_{ik} if x_i=x_k and -c_{ik} otherwise.
    I_star->push_back(k);
    d[k].second = true;
    for (auto iter = qi_.get_nonzero_begin(k); iter != qi_.get_nonzero_end(k);
	 ++iter) {
      if (assignments_[iter->first] == assignments_[k]) {
	d[iter->first].first += 2.0 * iter->second;
      } else {
	d[iter->first].first -= 2.0 * iter->second;
      }
    }
  }
}

void Palubeckis2004bSolution::SteepestAscent(const std::vector<int> I_star) {
  // Build a bitmap of set inclusion in I_star. Will be updated as we remove
  // elements.
  std::vector<bool> in_I_star(N_, false);
  for (int ct=0; ct < I_star.size(); ++ct) {
    in_I_star[I_star[ct]] = true;
  }

  // We will select some subset of indices in I_star to flip. They will be added
  // to this vector, and then we'll flip all of them at the end.
  std::vector<int> to_flip;

  // Step 1: Build h1_i and h2_i for each i in I_star
  std::vector<double> h1(N_, 0.0);
  std::vector<double> h2(N_, 0.0);
  for (int ct=0; ct < I_star.size(); ++ct) {
    int i = I_star[ct];
    h1[i] = diff_weights_[i];
    for (auto iter = qi_.get_nonzero_begin(i); iter != qi_.get_nonzero_end(i);
	 ++iter) {
      if (in_I_star[iter->first]) {
	// c'_ij is c_ij if x[i]=x[j] and is -c_ij if x[i] != x[j]
	if (assignments_[i] == assignments_[iter->first]) {
	  h2[i] += 2.0 * iter->second;
	} else {
	  h2[i] -= 2.0 * iter->second;
	}
      }
    }
  }

  // Step 2: Loop for |I_star|
  for (int iter_count=0; iter_count < I_star.size(); ++iter_count) {
    // Step 2.1: Initialize V1 and V2 (and a few other vars)
    double V1 = -std::numeric_limits<double>::max();
    double V2 = -std::numeric_limits<double>::max();
    int j = -1;
    int v = -1;

    // Step 2.2: Loop through the remaining elements if I_star
    for (int ct=0; ct < I_star.size(); ++ct) {
      int i = I_star[ct];
      if (!in_I_star[i]) {
	continue;  // This element has already been handled
      }
      
      // Step 2.2.1: Initialize q1 and q2
      double q1 = 2.0 * h1[i] + h2[i];
      double q2 = h1[i];

      // Step 2.2.2: Set the r value
      int r;
      if (q1 > 0.0 || (q1 == 0.0 && q2 >= 0.0)) {
	r = 1;
      } else {
	r = 0;
	q1 = -q1;
	q2 = -q2;
      }

      // Step 2.2.3: Update if best move found
      if (q1 > V1 || (q1 == V1 && q2 > V2)) {
	V1 = q1;
	V2 = q2;
	j = i;
	v = r;
      }
    }

    // Step 2.3: Update to_flip, I_star, h1, and h2 based on j
    if (v == 1) {
      to_flip.push_back(j);
    }
    in_I_star[j] = false;
    for (auto iter = qi_.get_nonzero_begin(j); iter != qi_.get_nonzero_end(j);
	 ++iter) {
      if (in_I_star[iter->first]) {
	// c'_ij is c_ij if x[i]=x[j] and is -c_ij if x[i] != x[j]
	if (assignments_[j] == assignments_[iter->first]) {
	  h2[iter->first] -= 2.0 * iter->second;
	  if (v == 1) {
	    h1[iter->first] += 2.0 * iter->second;
	  }
	} else {
	  h2[iter->first] += 2.0 * iter->second;
	  if (v == 1) {
	    h1[iter->first] -= 2.0 * iter->second;
	  }
	}
      }
    }
  }

  // Step 6 of MST2 (combined in here): flip all selected elements in I_star
  for (int ct=0; ct < to_flip.size(); ++ct) {
    UpdateCutValues(to_flip[ct]);
  }
}

Palubeckis2004bMST2::Palubeckis2004bMST2(const QUBOInstance& qi,
					 double runtime_limit, bool validation,
                                         QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z1;
  int Z2;
  if (qi.get_size() > 500) {
    Z1 = 25000;
    Z2 = 10000;
  } else {
    Z1 = 10000;
    Z2 = 2500;
  }
  int z1max = std::max(500000, qi.get_size() * Z1);
  int z2max = std::max(500000, qi.get_size() * Z2);
  double alpha = 0.4;
  int n_prime = std::max<int>(10, (int)floor(alpha * qi.get_size()));

  // Step 1: Generate a random solution, and set it as the best encountered.
  Palubeckis2004bSolution x = QUBOSolution::RandomSolution(qi, this);
  double best_objective = x.get_weight();

  // Step 2: Improve this initial random solution with STS. If termination
  // criterion is met, stop execution.
  x.STS(&best_objective, z1max);

  // Step 3: Loop while termination criterion not met
  while (Report(x)) {
    // Step 4: Run SELECT_VARIABLES to select indices in I_star
    std::vector<int> I_star;
    x.SelectVariables(n_prime, &I_star);

    // Steps 5+6: Run STEEPEST_ASCENT to flip some subset of the indices in
    // I_star, updating the solution variables based on these flipped values.
    x.SteepestAscent(I_star);

    // Step 7: Run STS on x, maintaining x as the current solution as opposed
    // to taking the best solution encountered.
    x.STS(&best_objective, z2max);
  }
}

Palubeckis2004bSolution::Palubeckis2004bSolution(const QUBOInstance& qi, int mu,
						 QUBOHeuristic *heuristic) :
  QUBOSolution(qi, heuristic) {
  // Step 1: Initialize g1 and g2 vectors
  std::vector<double> g1 = qi.get_lin();  // g1 initialized to linear terms
  std::vector<double> g2 = std::vector<double>(qi.get_size(), 0.0);
  for (auto iter = qi.get_all_nonzero_begin(); iter != qi.get_all_nonzero_end();
       ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double q_ij = iter->second;
    g2[i] += 2.0 * q_ij;
    g2[j] += 2.0 * q_ij;
  }
  // selected -- whether a variable has been removed from I
  std::vector<bool> selected = std::vector<bool>(qi.get_size(), false);

  for (int remaining = qi.get_size(); remaining > 0; --remaining) {
    // Step 2: If more than mu variables remain, select the mu of them with the
    // largest absolute value of 2g1+g2
    std::vector<std::pair<double, int> > sorted;
    for (int i=0; i < qi.get_size(); ++i) {
      if (!selected[i]) {
	sorted.push_back(std::pair<double, int>(fabs(2*g1[i] + g2[i]), i));
      }
    }
    if (remaining > mu) {
      std::sort(sorted.begin(), sorted.end(),
		std::greater<std::pair<double, int> >());
    }

    // Step 3: Use roulette wheel selection to pick the node to be added
    std::vector<double> scores;
    int numSelect = std::min<int>(mu, sorted.size());
    for (int count=0; count < numSelect; ++count) {
      scores.push_back(sorted[count].first);
    }
    int index = Random::RouletteWheel(scores);
    int j = sorted[index].second;
    
    // Step 4: Assign j to the proper set based on its q_j value, and update
    // g1 and g2 values accordingly.
    if (2*g1[j] + g2[j] >= 0) {
      assignments_[j] = 1;
    }
    selected[j] = true;
    for (auto iter = qi.get_nonzero_begin(j); iter != qi.get_nonzero_end(j);
	 ++iter) {
      g2[iter->first] -= 2.0 * iter->second;
      if (assignments_[j]) {
	g1[iter->first] += 2.0 * iter->second;
      }
    }
  }

  // All variables have been selected, so fill in diff_weights_ and weight_
  PopulateFromAssignments();
}

Palubeckis2004bMST3::Palubeckis2004bMST3(const QUBOInstance& qi,
					 double runtime_limit, bool validation,
                                         QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z;
  if (qi.get_size() >= 2500) {
    Z = 10700;
  } else if (qi.get_size() >= 1000) {
    Z = 10200;
  } else {
    Z = 2500;
  }
  int zmax = std::max(500000, qi.get_size() * Z);
  int mu = 50;

  // Step 1: Set the best ever objective value to 0 (score from all-0 solution)
  double best_objective = 0.0;

  while (1) {
    // Step 2: Build vector with GET_START_POINT
    Palubeckis2004bSolution solution =
      Palubeckis2004bSolution::GetStartPoint(qi, mu, this);

    // Step 3: Apply STS to the solution
    solution.STS(&best_objective, zmax);

    // If termination criterion is met, stop execution.
    if (!Report(solution)) {
      break;
    }
  }
}

// RandomFromElite, the MST4 Step 4 constructive procedure
Palubeckis2004bSolution::Palubeckis2004bSolution(const std::vector<Palubeckis2004bSolution>& B,
						 double pPrime) :
  QUBOSolution(B[0]) {
  // Count number of times each variable is used in B
  std::vector<int> counts(N_, 0);
  for (auto iter = B.begin(); iter != B.end(); ++iter) {
    for (int i=0; i < N_; ++i) {
      counts[i] += iter->assignments_[i];
    }
  }

  // Select the variables for this solution
  for (int i=0; i < N_; ++i) {
    if (counts[i] != B.size() && counts[i] != 0) {
      // Some elite solutions different on this index
      double zeta;
      if (Random::RandDouble() <= pPrime) {
	zeta = Random::RandDouble();
      } else {
	zeta = 0.5;
      }

      int newAssignment = 0;
      if (zeta <= ((double)counts[i]) / B.size()) {
	newAssignment = 1;
      }

      if (newAssignment != assignments_[i]) {
	UpdateCutValues(i);
      }
    }
  }
}

Palubeckis2004bMST4::Palubeckis2004bMST4(const QUBOInstance& qi,
					 double runtime_limit, bool validation,
                                         QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z;
  if (qi.get_size() >= 2500) {
    Z = 11500;
  } else if (qi.get_size() >= 1000) {
    Z = 9700;
  } else {
    Z = 3100;
  }
  int zmax = std::max(500000, qi.get_size() * Z);
  double Delta = 0.2;
  int mStar = 100;
  double pPrime = 1.0;
  int tStar = 100;

  // Do random restarts using the original termination criterion; because this
  // algorithm using global state that will eventually converge it's wasteful to
  // continue indefinitely instead of using random restarts.
  while (true) {
    // Build sets
    std::vector<Palubeckis2004bSolution> B;  // Elite solutions from tabu search
    std::vector<Palubeckis2004bSolution> S;  // Set used to avoid repeats

    // Step 1: Generate a random solution, and set it as the best encountered.
    Palubeckis2004bSolution initial = QUBOSolution::RandomSolution(qi, this);
    double best_objective = initial.get_weight();
    S.push_back(initial);

    // Step 2: Apply STS_MOD to this initial solution
    initial.STSMod(&best_objective, zmax, &B, mStar, Delta);

    // Report initial solution, and exit if out of time
    if (!Report(initial)) {
      return;
    }

    for (int t=0; t < tStar; ++t) {
      // Step 4: Construct solution from set B
      Palubeckis2004bSolution solution =
        Palubeckis2004bSolution::RandomFromElite(B, pPrime);
      
      // Step 5: Check if this solution matches a member of S
      bool matches = false;
      for (auto iter = S.begin(); iter != S.end(); ++iter) {
        if (solution == *iter) {
          matches = true;
          break;
        }
      }
      
      // Step 6: If the solution didn't match (I'm pretty sure that's what
      // they're going for here, though it's not completely clear from the
      // pseudocode), then add this solution to S and run STS_MOD, which will
      // generate a new set B.
      if (!matches) {
        S.push_back(solution);
        solution.STSMod(&best_objective, zmax, &B, mStar, Delta);
      }
      
      // Check termination criterion
      if (!Report(solution)) {
        return;
      }
    }
  }
}

Palubeckis2004bInstance::Palubeckis2004bInstance(const QUBOInstance& qi,
						 const Palubeckis2004bSolution& x) :
  QUBOInstance(qi) {
  // Parameters
  double psel = 0.4;
  int delta = 1;
  
  // Let's start by perturbing the non-zero values from all_nonzero_
  for (auto iter = get_mutable_all_nonzero_begin();
       iter != get_mutable_all_nonzero_end(); ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    if (Random::RandDouble() <= psel) {
      if (x.get_assignments()[i] == 1 && x.get_assignments()[j] == 1) {
	iter->second -= Random::RandInt(0, delta) / 2.0;
      } else {
	iter->second += Random::RandInt(0, delta) / 2.0;
      }
    }
  }

  // Now, let's refresh the nonzero_ vectors
  for (int ct=0; ct < nonzero_.size(); ++ct) {
    nonzero_[ct].clear();
  }
  for (auto iter = get_all_nonzero_begin(); iter != get_all_nonzero_end();
       ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double q_ij = iter->second;
    nonzero_[i].push_back(std::pair<int, double>(j, q_ij));
    nonzero_[j].push_back(std::pair<int, double>(i, q_ij));
  }
}

Palubeckis2004bSolution::Palubeckis2004bSolution(const QUBOInstance& perturbed,
						 const Palubeckis2004bSolution& x) :
  QUBOSolution(perturbed, x.heuristic_) {
  // Use the assignments from the passed solution, but assign diff_weights_ and
  // weight_ with the passed instance.
  assignments_ = x.assignments_;
  PopulateFromAssignments();
}

Palubeckis2004bSolution::Palubeckis2004bSolution(const Palubeckis2004bSolution& x,
						 int ignored) :
  QUBOSolution(x.qi_, x.heuristic_) {
  // Step 0: Assign parameters
  int Z3;
  if (qi_.get_size() >= 2500) {
    Z3 = 2250;
  } else if (qi_.get_size() >= 1000) {
    Z3 = 1750;
  } else {
    Z3 = 200;
  }
  int z3max = std::max(500000, qi_.get_size() * Z3);

  // Step 2: Construct a perturbed problem instance
  QUBOInstance perturbed = Palubeckis2004bInstance(qi_, x);

  // Step 3: Use the perturbed instance for this solution. We achieve this by
  // actually constructing a new object with the perturbed instance.
  Palubeckis2004bSolution y(perturbed, x);

  // Step 4: Run STS on new instance, tracking the best solution obtained in
  // variable best but not reporting new best solutions (because they have
  // the wrong objective due to the modified problem).
  Palubeckis2004bSolution best(y);
  double best_objective = best.weight_;
  y.STS(&best_objective, z3max, &best, false);

  // Step 5: Copy back assignments from best, and populate diff_weights_ and
  // weight_ from those assignment but the original QUBOInstance.
  assignments_ = best.assignments_;
  PopulateFromAssignments();
}

Palubeckis2004bMST5::Palubeckis2004bMST5(const QUBOInstance& qi,
					 double runtime_limit, bool validation,
                                         QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z1;
  int Z2;
  if (qi.get_size() > 500) {
    Z1 = 25000;
    Z2 = 10000;
  } else {
    Z1 = 10000;
    Z2 = 2100;
  }
  int z1max = std::max(500000, qi.get_size() * Z1);
  int z2max = std::max(500000, qi.get_size() * Z2);
  int tStar = 100;

  // Because this algorithm perturbs and then optimizes a single solution, we
  // will wrap it in random restarts.
  while (true) {
    // Step 1: Randomly generate solution, and set as best
    Palubeckis2004bSolution init = QUBOSolution::RandomSolution(qi, this);
    Palubeckis2004bSolution best(init);
    double best_objective = best.get_weight();
    
    // Step 2: Run STS on initial solution. best will contain the best solution
    // encountered during the course of the run, and init will contain the final
    // solution after the run (x* and x from the paper, respectively).
    init.STS(&best_objective, z1max, &best, true);
    
    // Step 3: Loop until termination criterion reached
    for (int t=0; t < tStar; ++t) {
      // Steps 4-5: Obtain a solution with PERTURB
      Palubeckis2004bSolution x = Palubeckis2004bSolution::Perturb(best);
      
      // Step 6: Apply STS, again updating the best-ever solution if improved
      x.STS(&best_objective, z2max, &best, true);
      
      // Check termination criterion
      if (!Report(best)) {
        return;
      }
    }
  }
}

Palubeckis2004bSTS::Palubeckis2004bSTS(const QUBOInstance& qi,
				       double runtime_limit, bool validation,
                                       QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Assign parameters (based on breakdown provided in results section)
  int Z;
  if (qi.get_size() >= 2500) {
    Z = 1650000;
  } else if (qi.get_size() >= 1000) {
    Z = 1400000;
  } else {
    Z = 380000;
  }
  int zmax = std::max(500000, qi.get_size() * Z);

  // Apply STS until termination, starting with random solution
  double best_objective = 0.0;
  while (true) {
    Palubeckis2004bSolution x = QUBOSolution::RandomSolution(qi, this);
    x.STS(&best_objective, zmax);
    if (!Report(x)) {
      break;
    }
  }
}
