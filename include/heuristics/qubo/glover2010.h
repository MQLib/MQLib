#ifndef HEURISTICS_QUBO_GLOVER_2010_H_
#define HEURISTICS_QUBO_GLOVER_2010_H_

#include <vector>
#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"
#include "util/random.h"

class Glover2010QUBOSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Glover2010QUBOSolution
  Glover2010QUBOSolution(const QUBOSolution &x);

  static Glover2010QUBOSolution
    PerturbationOperator(const Glover2010QUBOSolution& x,
			 const std::vector<int>& EliteFreq, int r,
			 const std::vector<int>& FlipFreq) {
    return Glover2010QUBOSolution(x, EliteFreq, r, FlipFreq);
  }

  // Run the tabu search procedure TS0, updating this solution. This function
  // populates the output argument FlipFreq with the frequency with which each
  // variable with flipped during the tabu search
  void TabuSearch(std::vector<int>* FlipFreq);

 private:
  // Run the perturbation operator (Sec. 2.3.2) from base solution x, with
  // provided elite frequencies and flip frequencies. r is the number of elite
  // solutions.
  Glover2010QUBOSolution(const Glover2010QUBOSolution& x,
			 const std::vector<int>& EliteFreq, int r,
			 const std::vector<int>& FlipFreq);
};

// Maintain a priority queue of elite solutions (with the worst solution at the
// top); it also maintains the freq of each variable being set in the elite set.
class Glover2010Elite {
 public:
  Glover2010Elite(const QUBOInstance& qi, int R);

  // Getters
  int size() const {  return EliteSol_.size();  }
  const Glover2010QUBOSolution& RandomSolution() const {
    return EliteSol_[Random::RandInt(0, EliteSol_.size()-1)];
  }
  const std::vector<int>& get_freq() const {  return EliteFreq_;  }

  // Add a solution if it isn't already in the elite set and either 1) The elite
  // set is not full size; or 2) It's better than the worst solution in the
  // elite set.
  void AddSolution(const Glover2010QUBOSolution& x);

 private:
  // Disable default constructor
  Glover2010Elite();

  // Maximum size of the elite set
  int R_;

  // Number of variables in the QUBOInstance being processed
  int N_;

  // EliteSol_ vector is maintained with make_heap, push_heap, and pop_heap;
  // this vector maintains the worst solution as the top of the heap.
  std::vector<Glover2010QUBOSolution> EliteSol_;

  // EliteFreq_ vector is the frequency of each variable being 1 in the 
  // EliteSol_ vector.
  std::vector<int> EliteFreq_;
};

class Glover2010 : public QUBOHeuristic {
 public:
  Glover2010(const QUBOInstance& qi, double runtime_limit, bool validation,
             QUBOCallback *qc);
};

#endif
