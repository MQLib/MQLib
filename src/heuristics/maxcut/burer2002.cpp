#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/maxcut/burer2002.h"
#include "problem/heuristic.h"
#include "util/random.h"

// Perform all first available 1-swaps. This is nearly identical to
// ExtendedSolution::AllFirst1Swap, except it has a tolerance.
void Burer2002Solution::All1Swap(double tolerance) {
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (int i=0; i < N_; ++i) {
      if (diff_weights_[i] > tolerance && ImprovingMove(i)) {
	UpdateCutValues(i);
	move_made = true;
	break;
      }
    }
  }
}

// Perform all first available 2-swaps. This is nearly identical to
// MaxCutSolution::AllFirst2Swap, except it has a tolerance.
void Burer2002Solution::All2Swap(double tolerance) {
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (auto iter = mi_.get_all_edges_begin();
	 iter != mi_.get_all_edges_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double w_ij = iter->second;
      double benefit = diff_weights_[i] + diff_weights_[j] -
	2.0 * assignments_[i] * assignments_[j] * w_ij;
      if (benefit > tolerance && ImprovingMove(benefit)) {
	UpdateCutValues(i);
	UpdateCutValues(j);
	move_made = true;
	break;
      }
    }
  }
}


// Given a new theta vector, compute the element-wise cosine and sine, as well
// as weighted sum of the cosines of the angle differences between paired
// nodes (returned through command-line arguments). Return objective.
double Burer2002Solution::LoadNewTheta(const std::vector<double>& theta,
				       std::vector<double>* cos_theta,
				       std::vector<double>* sin_theta,
				       std::vector<double>* dH) {
  for (int i=0; i < N_; ++i) {
    (*cos_theta)[i] = cos(theta[i]);
    (*sin_theta)[i] = sin(theta[i]);
    (*dH)[i] = 0.0;
  }
  double objective = 0.0;
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    double w_ij = iter->second;
    double scaled_cos_diff = w_ij *
      ((*cos_theta)[i] * (*cos_theta)[j] + (*sin_theta)[i] * (*sin_theta)[j]);
    objective += scaled_cos_diff;
    (*dH)[i] -= scaled_cos_diff;
    (*dH)[j] -= scaled_cos_diff;
  }
  return objective;
}

// Run Procedure-CUT to generate a new solution
Burer2002Solution::Burer2002Solution(const MaxCutInstance& mi,
				     double w1norm,
				     std::vector<double>* theta,
				     MaxCutHeuristic *heuristic) :
  MaxCutSolution(mi, heuristic) {
  // ***Parameters for construction method
  // Stopping value for relative change in f(theta) when optimizing f
  const double f_tolerance = 1e-4;
  // Stopping value for normalized value of g(theta) when optimizing f
  const double g_tolerance = 1e-4;
  // Maximum number of backtracks before termination
  const int maxback = 10;
  // Backtrack multiplier
  const double tau = 0.5;
  // The "beta" value in the Armijo line-search (required fraction of 1st order
  // reduction) -- see http://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf
  // slide 9
  const double gamma = 0.01;
  // The initial step size in the line search
  const double alpha_init = 1.0;
  // The maximum step size in the line search
  const double max_alpha = 4.0;
  // Maximum number of gradient descent steps when optimizing f(theta)
  const int max_opt_iter = 200;

  // Minimize f(theta) using gradient descent with backtracking Armijo
  // line-search.

  // *** Setup variables
  // The backtracking alpha value for each iteration's line search
  double bt_alpha = alpha_init;
  // negative of cosine differences incident to each node
  std::vector<double> dH(N_);
  // element-wise cosine of the theta vector
  std::vector<double> cos_theta(N_);
  // element-wise sine of the theta vector
  std::vector<double> sin_theta(N_);
  // f: the current objective value of nonlinear optimization problem
  double f = LoadNewTheta(*theta, &cos_theta, &sin_theta, &dH);
  for (int opt_iter = 0; opt_iter < max_opt_iter; ++opt_iter) {
    // Compute the current gradient
    std::vector<double> g(N_, 0.0);  // Gradient of f(theta)
    for (auto iter = mi_.get_all_edges_begin();
	 iter != mi_.get_all_edges_end(); ++iter) {
      int i = iter->first.first;
      int j = iter->first.second;
      double w_ij = iter->second;
      double scaled_sin_diff =
	w_ij * (sin_theta[j]*cos_theta[i] - cos_theta[j]*sin_theta[i]);
      g[i] += scaled_sin_diff;
      // sin(theta[i] - theta[j]) = -sin(theta[j] - theta[i])
      g[j] -= scaled_sin_diff;
    }
    
    // Stop the gradient descent if the gradient is too close to 0
    double norm_gradient = 0.0;
    for (int ct=0; ct < N_; ++ct) {
      norm_gradient += g[ct] * g[ct];
    }
    if (norm_gradient / w1norm < g_tolerance) {
      break;
    }
    
    // Determine the descent direction via scaling; we determined this behavior
    // by actually looking at the circut code as this is not mentioned in the
    // Burer2002 paper
    std::vector<double> desc(N_);  // Descent direction
    double g_times_desc = 0.0;
    double divisor = 1.0;
    for (int ct=0; ct < N_; ++ct) {
      divisor = std::max(divisor, dH[ct]);
    }
    for (int ct=0; ct < N_; ++ct) {
      desc[ct] = -g[ct] / divisor;
      g_times_desc += desc[ct] * g[ct];
    }
    
    // Use backtracking Armijo line-search to determine a good step size
    std::vector<double> new_theta(N_);
    int numback;
    double recent_f = -1.0;
    for (numback=1; numback <= maxback; ++numback) {
      // Compute the new theta with step size alpha
      for (int ct=0; ct < N_; ++ct) {
	new_theta[ct] = (*theta)[ct] + bt_alpha * desc[ct];
      }
      
      // Update the cosine and sine vectors, as well as dH and the
      // objective. Exit on Armijo condition.
      recent_f = LoadNewTheta(new_theta, &cos_theta, &sin_theta, &dH);
      if (recent_f <= f + gamma * bt_alpha * g_times_desc) {
	break;
      }
      
      // Exponential scaling on step size
      bt_alpha *= tau;
    }
    double f_prev = f;
    f = recent_f;
    *theta = new_theta;  // Copy over the current theta
    
    // Stop the gradient descent if there was not enough change in the
    // objective value
    double rel_change = fabs(f - f_prev) / (1.0 + fabs(f_prev));
    if (rel_change < f_tolerance) {
      break;
    }
    
    // Step size update lifted from circut code; not mentioned in paper...
    if (numback <= 2) {
      if (max_alpha < 2.0 * bt_alpha) {
	bt_alpha = max_alpha;
      } else {
	bt_alpha *= 2.0;
      }
    }
  }
  
  // Modulo the angles to be between 0 and 2*PI, and add to a vector of
  // index/angle pairs. Sort on angle.
  std::vector<std::pair<double, int> > angles;
  for (int ct=0; ct < N_; ++ct) {
    (*theta)[ct] -= 2 * M_PI * floor((*theta)[ct] / (2*M_PI));
    angles.push_back(std::pair<double, int>((*theta)[ct], ct));
  }
  angles.push_back(std::pair<double, int>(2*M_PI, N_));
  std::sort(angles.begin(), angles.end());

  // Determine initial set inclusion, and setup variables to be updated
  // later.
  for (int ct=0; ct < N_; ++ct) {
    assignments_[ct] = -1;
  }
  std::vector<std::pair<double, int> >::iterator first_it = angles.begin();
  std::vector<std::pair<double, int> >::iterator second_it =
    angles.begin();
  while (second_it->first <= M_PI) {
    assignments_[second_it->second] = 1;  // The first set has angles in [0, pi]
    ++second_it;
  }

  // Now, first_it points at the beginning, and second_it points at the first
  // node with theta > pi.

  // Fill in diff_weights_ and weight_ from assignments_, and copy them over to
  // the temporary solution we'll be updating in the search for a better
  // solution.
  PopulateFromAssignments();
  double curr_weight = weight_;
  std::vector<int> curr_assignments(assignments_);
  std::vector<double> curr_diff_weights(diff_weights_);
  
  // Compute the optimal cut by exhaustively searching through the possible cuts
  while (1) {
    // Update the cut angle
    int update_index;
    if (first_it->first <= second_it->first - M_PI) {
      update_index = first_it->second;  // We will exclude this guy.
      ++first_it;
    } else {
      update_index = second_it->second;  // We will include this guy
      ++second_it;
    }
    if (update_index == N_) {
      break;  // We have exhausted all the possible cuts
    }
    
    // See if this is a new best cut; if so, copy as the current solution
    UpdateCutValues(update_index, &curr_assignments, &curr_diff_weights,
		    &curr_weight);
    if (BaseSolution::ImprovesOver(curr_weight, weight_)) {
      weight_ = curr_weight;
      assignments_ = curr_assignments;
      diff_weights_ = curr_diff_weights;
    }
  }
}

Burer2002::Burer2002(const MaxCutInstance& mi, double runtime_limit,
		     bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Parameters
  // Number of permitted non-improving perturbations to optimal theta before
  // search is stopped. This was set to a few different values in the
  // computational results of burer2002 (0, 4, and 8 for torus set; 0 and 10
  // for G-set; and 10 and 50 for spin-glass dataset) so we'll use 50 since
  // it was the most common choice (and an intermediate value).
  const int N = 10;
  // Whether we perform greedy 1- and 2-moves after 
  const int local_search = 1;
  // Amount of improvement required to do a 1-move
  const double one_move_tolerance = 0.01;
  // Amount of improvement required to do a 2-move
  const double two_move_tolerance = 0.1;
  // Perturbation in the range [-perturbation*PI, perturbation*PI]
  const double perturbation = 0.2;

  // Compute 1-norm of vec(W), the vectorization of the weight matrix.
  double w1norm = 0.0;
  for (auto iter=mi.get_all_edges_begin(); iter != mi.get_all_edges_end();
       ++iter) {
    w1norm += 2.0 * fabs(iter->second);  // Count both directions of edge
  }

  for (int iter=0; ; ++iter) {  // Random restart until termination criterion
    // Generate random starting set of angles
    std::vector<double> theta(mi.get_size());
    for (int ct=0; ct < mi.get_size(); ++ct) {
      theta[ct] = Random::RandDouble() * 2 * M_PI;
    }

    // Algorithm 1 from Section 4
    double best_weight = -std::numeric_limits<double>::max();
    int k = 0;  // # of runs without an improvement in the best solution
    while (k <= N) {
      // Rank2Cut minimizes f(theta) and then uses Procedure-CUT to get the
      // best cut associated with this theta.
      Burer2002Solution x = Burer2002Solution::Rank2Cut(mi, w1norm, &theta,
							this);

      // Perform local searches (all 1- and 2-moves better than the tolerance)
      if (local_search) {
	x.All1Swap(one_move_tolerance);
	x.All2Swap(two_move_tolerance);
      }

      // Check termination criterion (runtime on non-validation runs; iteration
      // count on validation runs).
      if (!Report(x, iter)) {
	return;
      }
      
      // Update counter keeping track of iterations without improvement
      if (x.ImprovesOver(best_weight)) {
	best_weight = x.get_weight();
	k = 0;
      } else {
	++k;
      }

      // Perturb the angles associated with the current solution
      for (int ct=0; ct < mi.get_size(); ++ct) {
	theta[ct] = M_PI / 2.0 * (1.0 - x.get_assignments()[ct]) +
	  perturbation * (2 * M_PI * Random::RandDouble() - M_PI);
      }
    }
  }
}
