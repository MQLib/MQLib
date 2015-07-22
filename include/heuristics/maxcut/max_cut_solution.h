#ifndef HEURISTICS_MAXCUT_MAX_CUT_SOLUTION_H_
#define HEURISTICS_MAXCUT_MAX_CUT_SOLUTION_H_

#include <vector>
#include "heuristics/extended_solution.h"
#include "heuristics/maxcut/max_cut_partial_solution.h"
#include "problem/max_cut_heuristic.h"

class MaxCutSolution : public ExtendedSolution {
 public:
  // Copy over a MaxCutPartialSolution with no unassigned values
  MaxCutSolution(const MaxCutPartialSolution& x);

  // Initialize a completely random solution (p=0.5 for each vertex)
  static MaxCutSolution RandomSolution(const MaxCutInstance& mi,
				       MaxCutHeuristic *heuristic) {
    return MaxCutSolution(mi, heuristic, 0, 0);  // Private constructor
  }

  // Initialize a random solution (p for each vertex provided)
  static MaxCutSolution RandomSolution(const MaxCutInstance& mi,
				       const std::vector<double>& p,
				       MaxCutHeuristic *heuristic) {
    return MaxCutSolution(mi, p, heuristic);
  }

  // Initialize from vector of assignments
  MaxCutSolution(const std::vector<int>& assignments, const MaxCutInstance& mi,
                 MaxCutHeuristic *heuristic);

  // Perform all available 2-moves, at each step selecting the most valuable.
  void AllBest2Swap()  {  AllBest2Swap(0);  }

  // Perform all 2-moves, at each step selecting the first improving move
  void AllFirst2Swap()  {  AllFirst2Swap(0);  }

  // Assignment operator
  MaxCutSolution& operator=(const MaxCutSolution &rhs);

  // Copy constructor
  MaxCutSolution(const MaxCutSolution& x);

  // This function initializes diff_weights_ and weight_, given that
  // assignments_ is populated.
  void PopulateFromAssignments();

 protected:
  // Initialize mi and heuristic members but nothing else
  // ********************** IMPORTANT NOTE ***************************
  // diff_weights_ is not properly set after using this constructor, so you
  // must either set that vector manually or by using PopulateFromAssignments
  // *****************************************************************
  MaxCutSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
		 int initVal = -1);

  // Switch the set of update_index, updating, the assignments (x), the
  // diff_weights, and the objective.
  void UpdateCutValues(int update_index, std::vector<int>* x,
		       std::vector<double>* diff_weights,
		       double *objective) const;

  using ExtendedSolution::UpdateCutValues;  // Unhide single-argument version

  const MaxCutInstance& mi_;
  // The associated heuristic (for reporting purposes)
  MaxCutHeuristic *heuristic_;
  
  // Perform all available 2-moves, at each step selecting the most valuable.
  void AllBest2Swap(int startpos);

  // Perform all 2-moves, at each step selecting the first improving move
  void AllFirst2Swap(int startpos);

 private:
  // Initialize a completely random solution
  MaxCutSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
		 int ignored1, int ignored2);

  // Random solution given vertex probabilities
  MaxCutSolution(const MaxCutInstance& mi, const std::vector<double>& p,
		 MaxCutHeuristic *heuristic);
};

class FirstFixedMaxCutSolution : public MaxCutSolution {
 public:
  // Initialize a completely random solution (p=0.5 for each vertex except 1st)
  static FirstFixedMaxCutSolution RandomSolution(const MaxCutInstance& mi,
						 MaxCutHeuristic *heuristic,
						 int fixedVal) {
    return FirstFixedMaxCutSolution(mi, heuristic, fixedVal, 0);
  }

  // Initialize a random solution (p for each vertex provided)
  static FirstFixedMaxCutSolution RandomSolution(const MaxCutInstance& mi,
						 const std::vector<double>& p,
						 MaxCutHeuristic *heuristic,
						 int fixedVal) {
    return FirstFixedMaxCutSolution(mi, p, heuristic, fixedVal);
  }

  // Perform all available 1-moves, at each step selecting the most valuable.
  // We set local search start pos to 1 so we don't flip the first index.
  void AllBest1Swap() {  ExtendedSolution::AllBest1Swap(1);  }

  // Perform all 1-moves, at each step selecting the first improving move.
  // We set local search start pos to 1 so we don't flip the first index.
  void AllFirst1Swap() {  ExtendedSolution::AllFirst1Swap(1);  }

  // Perform all 1-moves, iteratively taking the first improving move on a
  // shuffled list of indices. We set local search start pos to 1 so we don't
  // flip the first index.
  void AllShuffle1Swap() {  ExtendedSolution::AllShuffle1Swap(1);  }

  // Perform all available 2-moves, at each step selecting the most valuable.
  // We set local search start pos to 1 so we don't flip the first index.
  void AllBest2Swap() {  AllBest2Swap(1);  }

  // Perform all 2-moves, at each step selecting the first improving move.
  // We set local search start pos to 1 so we don't flip the first index.
  void AllFirst2Swap() {  AllFirst2Swap(1);  }

  // Assignment operator
  FirstFixedMaxCutSolution& operator=(const FirstFixedMaxCutSolution &rhs);

  // Copy constructor
  FirstFixedMaxCutSolution(const FirstFixedMaxCutSolution& x);

 protected:
  // Initialize mi, heuristic, and fixedVal_ members but nothing else
  FirstFixedMaxCutSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
			   int fixedVal);

  // This function initializes diff_weights_ and weight_, given that
  // assignments_ is populated.
  void PopulateFromAssignments();

  // Switch the set of update_index, updating, the assignments (x), the
  // diff_weights, and the objective.
  void UpdateCutValues(int update_index, std::vector<int>* x,
		       std::vector<double>* diff_weights,
		       double *objective) const;

  using ExtendedSolution::UpdateCutValues;  // Unhide single-argument version
  using MaxCutSolution::AllBest2Swap;  // Unhide 1-argument version
  using MaxCutSolution::AllFirst2Swap;  // Unhide 1-argument version

  // fixedVal_ is the fixed value for the first node (1 or -1)
  int fixedVal_;

 private:
  // Initialize a completely random solution
  FirstFixedMaxCutSolution(const MaxCutInstance& mi, MaxCutHeuristic *heuristic,
			   int fixedVal, int ignored);

  // Random solution given vertex probabilities
  FirstFixedMaxCutSolution(const MaxCutInstance& mi,
			   const std::vector<double>& p,
			   MaxCutHeuristic *heuristic,
			   int fixedVal);
};

#endif
