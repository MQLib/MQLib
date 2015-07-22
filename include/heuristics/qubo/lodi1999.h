#ifndef HEURISTICS_QUBO_LODI_1999_H_
#define HEURISTICS_QUBO_LODI_1999_H_

#include <vector>
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class Lodi1999MinRange {
 public:
  // Initialize the globally fixed indices (CF_) and associated upper and lower
  // bounds (LB_ and UB_)
  Lodi1999MinRange(const QUBOInstance& qi);

  // Given a set of fixed indices, perform the MinRange function, updating to
  // indicate other indices that should be fixed.
  void MinRange(std::vector<int>* fixed, ExtendedSolution *x) const;

  // Return whether each variable is fixed globally (in CF)
  const std::vector<int>& get_CF() const {  return CF_;  }
  
 private:
  // Perform MinRange, taking as input the lower bounds and upper bounds for
  // each index and updating them based on the new indices that are fixed.
  void MinRangeInternal(std::vector<int>* fixed, std::vector<double>* LB,
                        std::vector<double>* UB, ExtendedSolution* x) const;

  // Update bounds LB and UB when fixing variable i to value val
  void FixVariable(int i, int val, std::vector<double>* LB,
                   std::vector<double>* UB) const;

  // Instance associated with this run
  const QUBOInstance& qi_;
  
  // Globally fixed indices and associated bounds for MinRange function
  std::vector<int> CF_;  // -1 for unfixed; 0/1 for fixed to that value
  std::vector<double> LB_;
  std::vector<double> UB_;
};

class Lodi1999Solution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Lodi1999
  Lodi1999Solution(const QUBOSolution &x);

  // Random values (except those fixed with CF) and local search
  static Lodi1999Solution RandomWithCF(const std::vector<int>& CF,
                                       const QUBOInstance& qi,
                                       QUBOHeuristic* heuristic) {
    return Lodi1999Solution(CF, qi, heuristic);
  }

  // Crossover operator, mutation, and local search
  static Lodi1999Solution Child(const Lodi1999Solution &father,
                                const Lodi1999Solution &mother,
                                const Lodi1999MinRange& mr) {
    return Lodi1999Solution(father, mother, mr);
  }

  // Perform local search
  void LS(const std::vector<int>& fixed);

 private:
  // Random
  Lodi1999Solution(const std::vector<int>& CF, const QUBOInstance& qi,
                   QUBOHeuristic* heuristic);

  // Crossover
  Lodi1999Solution(const Lodi1999Solution &father,
                   const Lodi1999Solution &mother,
                   const Lodi1999MinRange& mr);
};


class Lodi1999 : public QUBOHeuristic {
 public:
  Lodi1999(const QUBOInstance& qi, double runtime_limit, bool validation,
           QUBOCallback *qc);
};

#endif
