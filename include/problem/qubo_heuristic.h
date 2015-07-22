#ifndef PROBLEM_QUBO_HEURISTIC_H_
#define PROBLEM_QUBO_HEURISTIC_H_

#include <vector>
class QUBOSimpleSolution;
#include "heuristics/qubo/qubo_simple_solution.h"
#include "problem/heuristic.h"
#include "problem/qubo_instance.h"

class QUBOCallback {
 public:
  virtual bool Report(const QUBOSimpleSolution& solution, bool newBest,
                      double runtime) = 0;
  virtual bool Report(const QUBOSimpleSolution& solution, bool newBest,
                        double runtime, int iter) = 0;
};

class QUBOHeuristic : public Heuristic {
 public:
  QUBOHeuristic(const QUBOInstance& qi, double runtime_limit, bool validation,
                QUBOCallback *c);

  // Heuristic checks whether or not it should continue without reporting a
  // new solution or an iteration count.
  bool Report();

  // Heuristic checks whether or not it should continue without reporting
  // a new solution but reporting an iteration count.
  bool Report(int val);

  // Heuristic reports a solution, which may or may not be the best, and no
  // iteration count.
  bool Report(const BaseSolution& solution);

  // Heuristic reports a solution, which may or may not be the best, and
  // an iteration count.
  bool Report(const BaseSolution& solution, int iter);

  // Check if all best reported solutions had accurate objective values
  bool IsHistoryValid();

  // Returns const reference to the best solution encountered
  const QUBOSimpleSolution& get_best_solution() const {
    return past_solutions_[past_solutions_.size()-1];
  }

 protected:
  // Determine if the passed solution is a new best solution; if so, store it
  // appropriately and return true.
  bool NewBest(const BaseSolution& solution, double runtime);

  // Instance to be run
  const QUBOInstance& qi_;

  // If validation_ is false, we only maintain the current best solution in this
  // vector. Otherwise we maintain all best solutions.
  std::vector<QUBOSimpleSolution> past_solutions_;

  // Callback (NULL for no callback)
  QUBOCallback *qc_;
};

#endif
