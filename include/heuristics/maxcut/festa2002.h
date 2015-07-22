#ifndef HEURISTICS_MAXCUT_FESTA_2002_H_
#define HEURISTICS_MAXCUT_FESTA_2002_H_

#include "problem/max_cut_heuristic.h"
#include "heuristics/maxcut/max_cut_partial_solution.h"
#include "heuristics/maxcut/max_cut_solution.h"

class Festa2002PartialSolution : public MaxCutPartialSolution {
 public:
  // Initialize using the adaptive greedy function described in Sec. 2.1
  static
    Festa2002PartialSolution AdaptiveGreedySolution(const MaxCutInstance& mi,
						    double alpha,
						    const std::vector<std::pair<std::pair<int, int>, double> >& sorted,
						    MaxCutHeuristic *heuristic) {
    return Festa2002PartialSolution(mi, alpha, sorted, heuristic);
  }

  const std::vector<int>& S() const {  return S_; }
  const std::vector<int>& Sbar() const {  return Sbar_; }

 private:
  // Initialize using the adaptive greedy function described in Sec. 2.1
  Festa2002PartialSolution(const MaxCutInstance& mi, double alpha,
			   const std::vector<std::pair<std::pair<int, int>, double> >& sorted,
			   MaxCutHeuristic *heuristic);


  // Ordered vectors of elements in S and Sbar
  std::vector<int> S_;
  std::vector<int> Sbar_;
};

class Festa2002Solution : public MaxCutSolution {
 public:
  // Convert from MaxCutSolution to Festa2002Solution
  Festa2002Solution(const MaxCutSolution &x);

  // Convert from Festa2002PartialSolution to Festa2002Solution
  Festa2002Solution(const Festa2002PartialSolution &x);

  // Initialize using path relinking as described in Sec. 2.2
  static Festa2002Solution PathRelinkingSolution(const MaxCutInstance& mi,
						 const Festa2002Solution& x,
						 const Festa2002Solution& z,
						 MaxCutHeuristic *heuristic) {
    return Festa2002Solution(mi, x, z, heuristic);
  }

  // Initialize to a random node in the selected neighborhood.
  static Festa2002Solution NeighborhoodSolution(const MaxCutInstance& mi,
						const Festa2002Solution& base,
						int neighborhood,
						MaxCutHeuristic *heuristic) {
    return Festa2002Solution(mi, base, neighborhood, heuristic);
  }

  // Perform VNS with maximum neighborhood kMax.
  void VNS(int kMax);

  // Perform local search
  void LocalSearch();

 private:
  // Initialize using path relinking as described in Sec. 2.2
  Festa2002Solution(const MaxCutInstance& mi, const Festa2002Solution& x,
		    const Festa2002Solution& z, MaxCutHeuristic *heuristic);

  // Initialize to a random node in the selected neighborhood.
  Festa2002Solution(const MaxCutInstance& mi, const Festa2002Solution& base,
		    int neighborhood, MaxCutHeuristic *heuristic);

  // Ordered vectors of elements in S and Sbar
  std::vector<int> S_;
  std::vector<int> Sbar_;
};

class Festa2002Heuristic : public MaxCutHeuristic {
 public:
  Festa2002Heuristic(const MaxCutInstance& mi, double runtime_limit,
		     bool validation, MaxCutCallback *mc, bool grasp, bool pr,
                     bool vns);
};

class Festa2002G : public Festa2002Heuristic {
 public:
  Festa2002G(const MaxCutInstance &mi, double runtime_limit, bool validation,
             MaxCutCallback *mc);
};

class Festa2002GPR : public Festa2002Heuristic {
 public:
  Festa2002GPR(const MaxCutInstance &mi, double runtime_limit, bool validation,
               MaxCutCallback *mc);
};

class Festa2002VNS : public Festa2002Heuristic {
 public:
  Festa2002VNS(const MaxCutInstance &mi, double runtime_limit, bool validation,
               MaxCutCallback *mc);
};

class Festa2002VNSPR : public Festa2002Heuristic {
 public:
  Festa2002VNSPR(const MaxCutInstance &mi, double runtime_limit,
		 bool validation, MaxCutCallback *mc);
};

class Festa2002GVNS : public Festa2002Heuristic {
 public:
  Festa2002GVNS(const MaxCutInstance &mi, double runtime_limit, bool validation,
                MaxCutCallback *mc);
};

class Festa2002GVNSPR : public Festa2002Heuristic {
 public:
  Festa2002GVNSPR(const MaxCutInstance &mi, double runtime_limit,
		  bool validation, MaxCutCallback *mc);
};

#endif
