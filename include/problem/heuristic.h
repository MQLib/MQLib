#ifndef PROBLEM_HEURISTIC_H_
#define PROBLEM_HEURISTIC_H_

#include <sys/time.h>
#include <sstream>
#include <string>
#include <vector>

class Heuristic {
 public:
  Heuristic(double runtime_limit, bool validation);

  // Compute the runtime from start until now
  double Runtime();

  // Getters
  double get_best() const {  return best_;  }

  /* In this section we have various functions for checking if the heuristic
   *   should keep going or stop. In all cases, the following hierarchy is used
   *   to determine whether to keep going:
   *   1) If the user provided a callback, use that callback
   *   2) Runtime-based termination criterion
   */

  // Summary of reported solutions
  std::string History();
  virtual bool IsHistoryValid() = 0;

  // Heuristic is an abstract class, so it should have a virtual destructor
  virtual ~Heuristic() {};

 protected:
  // Is this a validation run? If so we will store all solutions and evaluate
  // them after the run is complete.
  bool validation_;

  struct timeval start_time_;
  double best_;
  double runtime_limit_;

  // Information about the history of new best solutions found by the heuristic.
  std::vector<double> past_solution_values_;
  std::vector<double> past_solution_times_;

 private:
  // Disable default constructor
  Heuristic();
};

#endif
