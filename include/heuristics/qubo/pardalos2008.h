#ifndef HEURISTICS_QUBO_PARDALOS_2008_H_
#define HEURISTICS_QUBO_PARDALOS_2008_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"
#include <vector>

class Pardalos2008QUBOSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Pardalos2008QUBOSolution
  Pardalos2008QUBOSolution(const QUBOSolution &x);

  // Generate a solution from a base solution x and set of generation probs
  static Pardalos2008QUBOSolution
    GenerateSolution(const Pardalos2008QUBOSolution x,
		     const std::vector<double>& probs, int k) {
    return Pardalos2008QUBOSolution(x, probs, k);
  }

  // Run the tabu search procedure (Figure 2, Sec. 2.3), returning multiple
  // locally optimal solutions through output parameter R.
  void TabuSearch(std::vector<Pardalos2008QUBOSolution>* R);

 private:
  // Run the generation procedure (Figure 3, Sec. 2.4) from base solution x
  // with generation probabilities probs.
  Pardalos2008QUBOSolution(const Pardalos2008QUBOSolution x,
			   const std::vector<double>& probs, int k);
};

// Maintain a priority queue of elite solutions previously encountered, also
// supporting an operation to remove all solutions sufficiently close to a
// supplied best solution.
class Pardalos2008Elite {
 public:
  Pardalos2008Elite(int Esize);

  // Getters
  int size() const {  return Elite_.size();  }
  void clear() {  Elite_.clear();  }
  const std::vector<Pardalos2008QUBOSolution>& get_solutions() const {
    return Elite_;
  }
  // What is the best elite solution?
  const Pardalos2008QUBOSolution& Best() const;

  // Add either one or multiple solutions to the elite set, adding solutions
  // only if they're in the top Esize_ best encountered.
  void AddSolution(const Pardalos2008QUBOSolution& x);
  void AddSolutions(const std::vector<Pardalos2008QUBOSolution>& slns);

  // Remove any elite solution that is within dp hamming distance of any
  // of the passed "best" solutions.
  void LimitByBests(const std::vector<Pardalos2008QUBOSolution>& bests);

 private:
  // Maximum size of the elite set
  int Esize_;

  // Elite_ vector is maintained with make_heap, push_heap, and pop_heap
  std::vector<Pardalos2008QUBOSolution> Elite_;
};

// Maintain the probabilities for each variable, enabling the efficient
// computation of ptilde_j. We maintain the numerator and denominator for
// E_{kj}^u.
class Pardalos2008Probs {
 public:
  // Initialize with an initial set of solutions
  Pardalos2008Probs(const std::vector<Pardalos2008QUBOSolution>& slns, int K,
		    const std::vector<double>& mu);

  // Add a set of solutions to the numerators and denominators of E_{kj}^u.
  void AddSolutions(const std::vector<Pardalos2008QUBOSolution>& slns);

  // Get the generation probabilities for each variable at temperature stage k;
  // store in output argument probs.
  void GetProbs(int k, std::vector<double>* probs) const;

 private:
  // There are K_+1 temperature stages; temps are in mu
  int K_;
  std::vector<double> mu_;

  // Number of variables in QUBO problem
  int N_;

  // Numerators and denominators of E_{kj}^u for u=0 and u=1. To access the
  // value for a given k, access at index k*N_ + j.
  std::vector<double> numerator0_;
  std::vector<double> numerator1_;
  std::vector<double> denominator0_;
  std::vector<double> denominator1_;

  // Frequencies of 1s and total number of solutions (needed to compute
  // p_j(0), which is the proportion of all j indexes with value 1).
  std::vector<int> freq1_;
  int freq_;
};

class Pardalos2008 : public QUBOHeuristic {
 public:
  Pardalos2008(const QUBOInstance& qi, double runtime_limit, bool validation,
               QUBOCallback *qc);
};

#endif
