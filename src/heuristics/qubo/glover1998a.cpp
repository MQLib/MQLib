#include <math.h>
#include <iostream>
#include <limits>
#include "heuristics/qubo/glover1998a.h"
#include "heuristics/qubo/qubo_solution.h"
#include "util/random.h"

Glover1998aSolution::Glover1998aSolution(const QUBOSolution &x) :
  QUBOSolution(x) {}

void Glover1998aSolution::doPhaseBegin(int matchVal, int k,
                                       const Glover1998aTabu& tabu,
                                       double PEN_R, double PEN_F) {
  // Convenience references
  const std::vector<int>& tabuR(tabu.get_tabuR());
  const std::vector<int>& tabuF(tabu.get_tabuF());

  // For the first k indices changed, use penalties from tabu values
  for (int iter=0; iter < k; ++iter) {
    double best_val = -std::numeric_limits<double>::max();
    int best_loc = -1;
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] == matchVal) {
        double improvement = diff_weights_[i] - PEN_R*tabuR[i] - PEN_F*tabuF[i];
        if (improvement > best_val) {
          best_val = improvement;
          best_loc = i;
        }
      }
    }
    if (best_loc < 0) {
      break;  // No more indices with value of matchVal
    }
    UpdateCutValues(best_loc);
  }

  // Use greedy best approach (only flipping if something equals matchVal)
  while (true) {
    double best_val = 0.0;
    int best_loc = -1;
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] == matchVal && diff_weights_[i] > best_val &&
          ImprovingMove(i)) {
        best_val = diff_weights_[i];
        best_loc = i;
      }
    }
    if (best_loc < 0) {
      break;  // No improving moves with value of matchVal
    }
    UpdateCutValues(best_loc);
  }
}

void Glover1998aSolution::doPhaseEnd(int matchVal, int span) {
  // Make the best available move for "span" iterations, limiting to moves
  // where the current value is matchVal
  for (int iter=0; iter < span; ++iter) {
    double best_val = 0.0;
    int best_loc = -1;
    for (int i=0; i < N_; ++i) {
      if (assignments_[i] == matchVal && diff_weights_[i] > best_val &&
          ImprovingMove(i)) {
        best_val = diff_weights_[i];
        best_loc = i;
      }
    }
    if (best_loc < 0) {
      break;  // No improving moves with value of matchVal
    }
    UpdateCutValues(best_loc);
  }
}

Glover1998aTabu::Glover1998aTabu(const QUBOInstance& qi) :
  t_(3),  // On page 337 they say this is the value used (could vary to 12...)
  N_(qi.get_size()),
  tabuR_(qi.get_size(), 0.0),
  tabuF_(qi.get_size(), 0.0),
  recent_(qi.get_size()*t_, 0.0),
  recent_pos_(0) {}

void Glover1998aTabu::CriticalEvent(const QUBOSolution& x) {
  // Subtract the element stored in recent_ at position recent_pos_ from tabuR_
  // Replace current recent_pos_ value in recent_ with x and update recent_pos_.
  // Add x to tabuR_ and tabuF_.
  int offset = recent_pos_*N_;
  const std::vector<int>& assignments = x.get_assignments();
  for (int i=0; i < N_; ++i) {
    tabuR_[i] += assignments[i] - recent_[offset+i];
    tabuF_[i] += assignments[i];
    recent_[offset+i] = assignments[i];
  }
  recent_pos_ = (recent_pos_+1) % t_;
}

void Glover1998a::transferPhase() {
  // Parameters
  const int p1 = 3;
  const int p2 = 7;
  const int intermediate_span_rate = 4;  // Hard-coded as 4 in section 4.1

  ++iter_span_;
  if (span_dir_increase_) {
    // Increasing span
    if (span_ <= p1) {
      if (iter_span_ > p2*span_) {
        ++span_;
        iter_span_ = 0;
      }
    } else {
      if (iter_span_ > intermediate_span_rate) {
        ++span_;
        iter_span_ = 0;
      }
    }
    if (span_ > p2) {
      span_ = p2;
      span_dir_increase_ = false;
    }
  } else {
    // Decreasing span
    if (span_ >= p1+1 && span_ <= p2) {
      if (iter_span_ > intermediate_span_rate) {
        --span_;
        iter_span_ = 0;
      }
    } else {
      if (iter_span_ > p2*span_) {
        --span_;
        iter_span_ = 0;
      }
    }
    if (span_ == 0) {
      span_ = 1;
      span_dir_increase_ = true;
      ++num_span_cycle_;
    }
  }
}

Glover1998a::Glover1998a(const QUBOInstance& qi, double runtime_limit,
			 bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc),
  tabu_(qi),
  N_(qi.get_size())
{
  // Parameters
  int cycle_limit = 5;  // Used for 50- and 100-node sets in paper
  if (N_ >= 500) {
    cycle_limit = 20;  // Used for 500-node sets in paper
  } else if (N_ >= 200) {
    cycle_limit = 10;  // Used for 200-node sets in paper
  }
  const int KMAX = 3;  // Top of page 339 suggests 2-4, so we pick 3
  int k_update_time = 2*3;  // This is thrown out as an example at the top of 339
  int max_span = 5;  // Maximum number of span cycles before restart

  // PEN_R is the maximum element of the A matrix and PEN_F=PEN_R/1000*#iter
  double PEN_R = 0.0;
  for (auto iter = qi.get_all_nonzero_begin(); iter != qi.get_all_nonzero_end();
       ++iter) {
    PEN_R = std::max<double>(PEN_R, iter->second);
  }
  for (int i=0; i < N_; ++i) {
    PEN_R = std::max<double>(PEN_R, qi.get_lin()[i]);
  }
  double PEN_F = PEN_R / 1000.0 * cycle_limit;

  while (true) {
    // Initialize to a random solution (though the paper uses an all-zero
    // solution, we'll use random instead to enable random restarts up to our
    // specified runtime limit).
    Glover1998aSolution x(QUBOSolution::RandomSolution(qi, this));
    
    // Initialize population parameters
    k_ = 1;
    iter_span_ = 0;
    span_ = 1;
    span_dir_increase_ = true;
    num_span_cycle_ = 0;
    
    // Loop until we've used our limit on the number of span cycles.
    for (int iter=0; num_span_cycle_ < max_span; ++iter) {
      // Loop through phases; first constructive (0->1) then destructive (1->0)
      for (int constructive=1; constructive >= 0; --constructive) {
        // Do phase up to critical event (penalties and matching vals differ
        // between constructive and destructive phases)
        x.doPhaseBegin(1-constructive, k_, tabu_, (2*constructive-1) * PEN_R,
                       (2*constructive-1) * PEN_F);
        
        // Report solution; exit if termination criterion. Report x to tabu list.
        if (!Report(x)) {
          return;  // Termination criterion met
        }
        tabu_.CriticalEvent(x);
        
        // Phase after critical event
        x.doPhaseEnd(1-constructive, span_);
        
        // Transfer phase
        transferPhase();
      }
      
      // Update k, the number of steps subject to the tabu penalties
      if ((iter+1) % k_update_time == 0) {
        k_ = (k_ == KMAX) ? 1 : (k_+1);
      }
    }
  }
}
