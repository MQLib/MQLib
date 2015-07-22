#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include "heuristics/maxcut/festa2002.h"
#include "problem/heuristic.h"
#include "util/random.h"

Festa2002Solution::Festa2002Solution(const MaxCutSolution &x) :
  MaxCutSolution(x) {
  // Fill in S_ and Sbar_ (they will be in sorted order)
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] == 1) {
      S_.push_back(i);
    } else {
      Sbar_.push_back(i);
    }
  }
}

// Convert from Festa2002PartialSolution to Festa2002Solution
Festa2002Solution::Festa2002Solution(const Festa2002PartialSolution &x) :
  MaxCutSolution(x),
  S_(x.S()),
  Sbar_(x.Sbar()) {}

Festa2002PartialSolution::Festa2002PartialSolution(const MaxCutInstance& mi,
						   double alpha,
						   const std::vector<std::pair<std::pair<int, int>, double> >& sorted,
						   MaxCutHeuristic *heuristic) :
  MaxCutPartialSolution(mi, heuristic) {
  // Though it's not clear from Section 2.1, build.f shows that the adaptive
  // greedy construction procedure starts by selecting an edge from a restricted
  // version of the sorted list of edges, restricting using alpha.
  int nrcl = std::min<int>(sorted.size(),
			   1 + (int)(sorted.size() * (1.0 - alpha)));
  int startEdge = Random::RandInt(0, nrcl - 1);
  UpdateCutValues(sorted[startEdge].first.first, 1);
  S_.push_back(sorted[startEdge].first.first);
  UpdateCutValues(sorted[startEdge].first.second, -1);
  Sbar_.push_back(sorted[startEdge].first.second);

  // Now, assign the remaining nodes.
  while (num_unassigned_ > 0) {
    // Determine the limits (w_min and w_max from paper)
    double this_best = -std::numeric_limits<double>::max();
    double this_worst = std::numeric_limits<double>::max();
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] != 0) {
	continue;  // i is already assigned
      }
      this_worst = std::min<double>(this_worst, gainS_[i]);
      this_worst = std::min<double>(this_worst, gainNS_[i]);
      this_best = std::max<double>(this_best, gainS_[i]);
      this_best = std::max<double>(this_best, gainNS_[i]);
    }

    // Determine cutoff (mu in paper) and find all uninserted nodes meeting
    // cutoff. The bool value indicates if they're to be added to S.
    double cutoff = this_worst + alpha * (this_best - this_worst) - 0.000001;
    std::vector<std::pair<int, bool> > made_cutoff;
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] == 0 && gainS_[i] >= cutoff) {
	made_cutoff.push_back(std::pair<int, bool>(i, true));
      }
      if (assignments_[i] == 0 && gainNS_[i] >= cutoff) {
	made_cutoff.push_back(std::pair<int, bool>(i, false));
      }
    }

    // Randomly select and add a valid node, updating sigma_S and sigma_NS
    int selectId = Random::RandInt(0, made_cutoff.size()-1);
    if (made_cutoff[selectId].second) {
      UpdateCutValues(made_cutoff[selectId].first, 1);
      S_.push_back(made_cutoff[selectId].first);
    } else {
      UpdateCutValues(made_cutoff[selectId].first, -1);
      Sbar_.push_back(made_cutoff[selectId].first);
    }
  }
}

Festa2002Solution::Festa2002Solution(const MaxCutInstance& mi,
				     const Festa2002Solution& x,
				     const Festa2002Solution& z,
				     MaxCutHeuristic *heuristic) :
  MaxCutSolution(x),
  S_(x.S_),
  Sbar_(x.Sbar_) {
  
  // Compute indices of symmetric difference between two linking solutions
  std::vector<int> diff;  // nodes that could switch sides
  int num_different = x.SymmetricDifference(z, &diff);

  // y is the current linking solution, and the class members contain the
  // best encountered to date. We start from solution x.
  std::vector<int> y_assignments = assignments_;
  std::vector<double> y_diff_weights = diff_weights_;
  double y_weight = x.get_weight();

  bool updated = false;
  for (int num_rem=num_different; num_rem >= 1; --num_rem) {
    // Determine the most valuable move in the current y solution (from
    // among solutions in the symmetric difference set).
    double best_diff = y_diff_weights[diff[0]];
    int best_index = 0;
    for (int i=1; i < diff.size(); ++i) {
      if (y_diff_weights[diff[i]] > best_diff) {
	best_diff = y_diff_weights[diff[i]];
	best_index = i;
      }
    }

    // Update y, flipping the selected vertex
    UpdateCutValues(diff[best_index], &y_assignments, &y_diff_weights,
		    &y_weight);

    // Remove the replaced index from the vector of differences
    diff[best_index] = diff[diff.size() - 1];
    diff.resize(diff.size() - 1);

    // If y is the best encountered middle solution, copy it over
    if (BaseSolution::ImprovesOver(y_weight, weight_)) {
      weight_ = y_weight;
      assignments_ = y_assignments;
      diff_weights_ = y_diff_weights;
      updated = true;
    }
  }

  // If we updated the solution, update S_ and Sbar_, as well
  if (updated) {
    std::vector<int> newS;
    std::vector<int> extraS;
    std::vector<int> newSbar;
    std::vector<int> extraSbar;
    for (int i=0; i < S_.size(); ++i) {
      if (assignments_[S_[i]] == 1) {
	newS.push_back(S_[i]);
      } else {
	extraSbar.push_back(S_[i]);
      }
    }
    for (int i=0; i < Sbar_.size(); ++i) {
      if (assignments_[Sbar_[i]] == -1) {
	newSbar.push_back(Sbar_[i]);
      } else {
	extraS.push_back(Sbar_[i]);
      }
    }
    std::copy(extraS.begin(), extraS.end(), std::back_inserter(newS));
    std::copy(extraSbar.begin(), extraSbar.end(), std::back_inserter(newSbar));
    S_ = newS;
    Sbar_ = newSbar;
  }
}

Festa2002Solution::Festa2002Solution(const MaxCutInstance& mi,
				     const Festa2002Solution& base,
				     int neighborhood,
				     MaxCutHeuristic *heuristic) :
  MaxCutSolution(base),
  S_(base.S_),
  Sbar_(base.Sbar_) {

  // Select k vertices to flip
  std::vector<int> to_flip(mi.get_size(), 0);
  for (int i=0; i < neighborhood; ++i) {
    int index;
    while (1) {
      index = Random::RandInt(0, mi.get_size()-1);
      if (!to_flip[index]) {
	to_flip[index] = 1;
	break;  // Found a new node to flip
      }
    }
  }

  // Adjust S_ and Sbar_ based on elements in to_flip -- first go through to
  // copy non-flipped values and then go through to add in flipped values at
  // end of respective vectors.
  std::vector<int> newS;
  std::vector<int> newSbar;
  for (int i=0; i < S_.size(); ++i) {
    if (!to_flip[S_[i]]) {
      newS.push_back(S_[i]);
    }
  }
  for (int i=0; i < Sbar_.size(); ++i) {
    if (!to_flip[Sbar_[i]]) {
      newSbar.push_back(Sbar_[i]);
    }
  }
  for (int i=0; i < N_; ++i) {
    if (to_flip[i]) {
      if (assignments_[i]) {
	newSbar.push_back(i);
      } else {
	newS.push_back(i);
      }
    }
  }
  S_ = newS;
  Sbar_ = newSbar;

  // Actually flip the vertices
  for (int i=0; i < mi.get_size(); ++i) {
    if (to_flip[i]) {
      UpdateCutValues(i);
    }
  }
}

void Festa2002Solution::VNS(int kMax) {
  int k = 1;
  int iter = 0;
  while (k <= kMax) {
    ++iter;
    if (iter % 1000 == 0) {
      // Terminate if the runtime limit has been reached, but report the iteration
      // count to be -1 so we won't terminate on a validation run.
      if (!heuristic_->Report(*this, -1)) {
	break;  // Out of time
      }
    }
    Festa2002Solution neighbor =
      Festa2002Solution::NeighborhoodSolution(mi_, *this, k, heuristic_);
    neighbor.LocalSearch();
    if (neighbor.ImprovesOver(*this)) {
      assignments_ = neighbor.assignments_;
      diff_weights_ = neighbor.diff_weights_;
      weight_ = neighbor.weight_;
      k = 1;
    } else {
      ++k;
    }
  }
}

void Festa2002Solution::LocalSearch() {
  // Maintain an ordered list of all elements in S and Sbar and go through
  // them in order, checking for improving moves
  if (S_.size() == 0 && Sbar_.size() == 0) {
    std::cout << "Missing S_ and Sbar_" << std::endl;
    exit(0);
  }

  std::vector<int> tmp;

  bool move_made = true;
  while (move_made) {
    move_made = false;

    // Loop through ordered lists, checking one from each list (this is not
    // well described in the paper, but is based on distributed code local.f)
    int Spos = 0;
    int Sbarpos = 0;
    bool moveInS = false;
    while (Spos < S_.size() || Sbarpos < Sbar_.size()) {
      if (Spos < S_.size()) {
	int var = S_[Spos];
	++Spos;
	if (ImprovingMove(var)) {
	  UpdateCutValues(var);
	  move_made = true;
	  moveInS = true;
	  break;
	}
      }

      if (Sbarpos < Sbar_.size()) {
	int var = Sbar_[Sbarpos];
	++Sbarpos;
	if (ImprovingMove(var)) {
	  UpdateCutValues(var);
	  move_made = true;
	  break;
	}
      }
    }

    // If a move was made, update the ordered lists accordingly
    if (move_made) {
      if (moveInS) {
	// Update Sbar_, adding moved element to back
	if (Sbarpos < Sbar_.size()) {
	  tmp.clear();
	  std::copy(Sbar_.begin() + Sbarpos, Sbar_.end(),
		    std::back_inserter(tmp));
	  std::copy(Sbar_.begin(), Sbar_.begin() + Sbarpos,
		    std::back_inserter(tmp));
	  tmp.push_back(S_[Spos-1]);
	  Sbar_ = tmp;
	} else {
	  Sbar_.push_back(S_[Spos-1]);
	}

	// Update S_, removing moved element
	if (Spos < S_.size()) {
	  tmp.clear();
	  std::copy(S_.begin() + Spos, S_.end(), std::back_inserter(tmp));
	  if (Spos >= 2) {
	    std::copy(S_.begin(), S_.begin() + Spos - 1,
		      std::back_inserter(tmp));
	  }
	  S_ = tmp;
	} else {
	  S_.pop_back();
	}
      } else {
	// Update S_, adding moved element to back
	if (Spos < S_.size()) {
	  tmp.clear();
	  std::copy(S_.begin() + Spos, S_.end(), std::back_inserter(tmp));
	  std::copy(S_.begin(), S_.begin() + Spos, std::back_inserter(tmp));
	  tmp.push_back(Sbar_[Sbarpos-1]);
	  S_ = tmp;
	} else {
	  S_.push_back(Sbar_[Sbarpos-1]);
	}

	// Update Sbar_, removing moved element
	if (Sbarpos < Sbar_.size()) {
	  tmp.clear();
	  std::copy(Sbar_.begin() + Sbarpos, Sbar_.end(),
		    std::back_inserter(tmp));
	  if (Sbarpos >= 2) {
	    std::copy(Sbar_.begin(), Sbar_.begin() + Sbarpos - 1,
		      std::back_inserter(tmp));
	  }
	  Sbar_ = tmp;
	} else {
	  Sbar_.pop_back();
	}
      }
    }
  }
}

Festa2002G::Festa2002G(const MaxCutInstance& mi, double runtime_limit,
		       bool validation, MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, true, false, false) {}

Festa2002GPR::Festa2002GPR(const MaxCutInstance& mi, double runtime_limit,
                           bool validation, MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, true, true, false) {}

Festa2002VNS::Festa2002VNS(const MaxCutInstance& mi, double runtime_limit,
                           bool validation, MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, false, false, true) {}

Festa2002VNSPR::Festa2002VNSPR(const MaxCutInstance& mi, double runtime_limit,
                               bool validation, MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, false, true, true) {}

Festa2002GVNS::Festa2002GVNS(const MaxCutInstance& mi,
			     double runtime_limit, bool validation,
                             MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, true, false, true) {}

Festa2002GVNSPR::Festa2002GVNSPR(const MaxCutInstance& mi, double runtime_limit,
                                 bool validation, MaxCutCallback* mc) :
  Festa2002Heuristic(mi, runtime_limit, validation, mc, true, true, true) {}

Festa2002Heuristic::Festa2002Heuristic(const MaxCutInstance& mi,
				       double runtime_limit, bool validation,
                                       MaxCutCallback* mc, bool grasp, bool pr,
                                       bool vns) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Parameters
  int MaxElite = 30;
  // The minimum symmetric difference to be considered "different enough"
  // from a member of the elite set.
  int sufficiently_different = (int)ceil(0.01 * mi.get_size());
  int kMax = 15;
  if (vns && !grasp) {
    kMax = 100;
  }
  if (kMax > mi.get_size()) {
    kMax = mi.get_size();
  }

  // If we're using GRASP, get a list of edges, sorted by weight
  std::vector<std::pair<std::pair<int, int>, double> > sorted;
  if (grasp) {
    mi.GetSortedEdges(&sorted);
  }

  std::vector<Festa2002Solution> elite;
  int first = 1;
  for (int i=0; ; ++i) {  // Iterate until termination criterion breaks loop
    Festa2002Solution solution = grasp ?
      Festa2002Solution(Festa2002PartialSolution::AdaptiveGreedySolution(mi, Random::RandDouble(), sorted, this)) :
      Festa2002Solution(MaxCutSolution::RandomSolution(mi, this));
    if (!vns) {
      solution.LocalSearch();
    } else {
      solution.VNS(kMax);
    }
    // If there is path relinking, do it now.
    if (pr) {
      if (first) {
	elite.push_back(solution);
        // Terminate if runtime limit is reached (non-validation run) or the
        // iteration limit has been reached (validation run).
	if (!Report(solution, i)) {
	  break;
	}
      } else {
	int idx = Random::RandInt(0, elite.size() - 1);
	Festa2002Solution pr_solution =
	  Festa2002Solution::PathRelinkingSolution(mi, solution, elite[idx],
						   this);
        // Terminate if runtime limit is reached (non-validation run) or the
        // iteration limit has been reached (validation run).
	if (!Report(pr_solution, i)) {
	  break;
	}
	if (elite.size() < MaxElite) {
	  elite.push_back(pr_solution);
	} else {
	  int too_close = 0;  // Too close to one of the elite solutions
	  int beaten = 0;  // Whether solution is beaten by an elite solution
	  int winner = 0;  // Whether solution beats an elite solution
	  double worst_obj = elite[0].get_weight();
	  int worst_index = 0;
	  for (int i=0; i < MaxElite; ++i) {
	    if (elite[i].get_weight() < worst_obj) {
	      worst_obj = elite[i].get_weight();
	      worst_index = i;
	    }
	    if (elite[i].ImprovesOver(pr_solution)) {
	      beaten = 1;
	    } else if (pr_solution.ImprovesOver(elite[i])) {
	      winner = 1;
	    }
	    if (pr_solution.SymmetricDifference(elite[i]) <
		sufficiently_different) {
	      too_close = 1;
	    }
	    if (too_close && beaten) {
	      break;  // Can't be selected as a new elite solution.
	    }
	  }
	  if (!beaten || (winner && !too_close)) {
	    elite[worst_index] = pr_solution;
	  }
	}
      }
    } else {
      // If no path relinking is used, just check if we have a best solution.
      // Terminate if runtime limit is reached (non-validation run) or the
      // iteration limit has been reached (validation run).
      if (!Report(solution, i)) {
	break;  // Termination criterion
      }
    }
    first = 0;
  }
}
