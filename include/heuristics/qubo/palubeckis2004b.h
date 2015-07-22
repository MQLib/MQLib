#ifndef HEURISTICS_QUBO_PALUBECKIS_2004B_H_
#define HEURISTICS_QUBO_PALUBECKIS_2004B_H_

#include "heuristics/qubo/qubo_solution.h"
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"
#include <vector>

// NOTE: In the paper Palubeckis2004b, for a given solution x the authors
// describe a procedure in which the problem instance is adjusted (via a simple
// variable substitution) to one for which that solution is all 0s. This enables
// them to ignore the interaction terms between variables and only focus on
// the (adjusted) linear terms. After this adjustment, the linear terms are
// exactly equal to the diff_weights_ values maintained by the solution class.
// Therefore, instead of adjusting the problem instances, we just maintain the
// diff_weights_ as usual.

class Palubeckis2004bSolution : public QUBOSolution {
 public:
  // Convert from QUBOSolution to Palubeckis2004bSolution
  Palubeckis2004bSolution(const QUBOSolution &x);

  // Runs construction procedure GET_START_POINT, used in MST3
  static Palubeckis2004bSolution GetStartPoint(const QUBOInstance& qi, int mu,
					       QUBOHeuristic *heuristic) {
    return Palubeckis2004bSolution(qi, mu, heuristic);
  }

  // Runs construction procedure from Step 4 of MST4
  static Palubeckis2004bSolution
    RandomFromElite(const std::vector<Palubeckis2004bSolution>& B,
		    double pPrime) {
    return Palubeckis2004bSolution(B, pPrime);
  }

  // Runs the PERTURB procedure, used in MST5
  static Palubeckis2004bSolution Perturb(const Palubeckis2004bSolution& x) {
    return Palubeckis2004bSolution(x, 0);
  }

  // Runs the SELECT_VARIABLES procedure, used in MST2. This generates n_prime
  // indices, which are stored in the output variable I_star.
  void SelectVariables(int n_prime, std::vector<int>* I_star) const;

  // Runs the STEEPEST_ASCENT procedure, used in MST2, as well as Step 6 of
  // MST2, which flips the selected elements of I_star.
  void SteepestAscent(const std::vector<int> I_star);

  // best_objective is a pointer to the best objective value ever found in the
  // master problem, Z, mStar, and Delta are parameters, and B is a vector to be
  // filled with elite solutions found during the tabu search (pass null for B
  // to run a standard STS method). best is enables the updates of the best
  // ever solution (NULL for no updates), and reportBest is true if you want
  // to report to the master heuristic and false if not.
  void STSMod(double *best_objective, int zmax,
	      std::vector<Palubeckis2004bSolution>* B, int mStar, double Delta,
	      Palubeckis2004bSolution* best = NULL, bool reportBest = true);

  // best_objective is a pointer to the best objective value ever found in the
  // master problem, and Z is a parameter.
  void STS(double* best_objective, int zmax,
	   Palubeckis2004bSolution *best = NULL, bool reportBest = true)  {
    STSMod(best_objective, zmax, NULL, 0, 0.0, best, reportBest);
  }

  // z is a pointer to the diff_weight_ check counter. LocalSearch performs
  // AllFirst1Swap updates.
  void LocalSearch(int *z);

 private:
  Palubeckis2004bSolution(const QUBOInstance& qi, int mu,
			  QUBOHeuristic *heuristic);

  Palubeckis2004bSolution(const std::vector<Palubeckis2004bSolution>& B,
			  double pPrime);

  Palubeckis2004bSolution(const Palubeckis2004bSolution& x, int ignored);

  // Create a solution with the same assignments as the passed solution, but
  // the weights from the passed QUBOInstance.
  Palubeckis2004bSolution(const QUBOInstance& perturbed,
			  const Palubeckis2004bSolution& x);
};

// Perturbs the QUBOInstance as in Step 2 of PERTURB() from MST5.
class Palubeckis2004bInstance : public QUBOInstance {
 public:
  Palubeckis2004bInstance(const QUBOInstance& qi,
			  const Palubeckis2004bSolution& x);
};

class Palubeckis2004bMST1 : public QUBOHeuristic {
 public:
  Palubeckis2004bMST1(const QUBOInstance& qi, double runtime_limit,
		      bool validation, QUBOCallback *qc);
};

class Palubeckis2004bMST2 : public QUBOHeuristic {
 public:
  Palubeckis2004bMST2(const QUBOInstance& qi, double runtime_limit,
		      bool validation, QUBOCallback *qc);
};

class Palubeckis2004bMST3 : public QUBOHeuristic {
 public:
  Palubeckis2004bMST3(const QUBOInstance& qi, double runtime_limit,
		      bool validation, QUBOCallback *qc);
};

class Palubeckis2004bMST4 : public QUBOHeuristic {
 public:
  Palubeckis2004bMST4(const QUBOInstance& qi, double runtime_limit,
		      bool validation, QUBOCallback *qc);
};

class Palubeckis2004bMST5 : public QUBOHeuristic {
 public:
  Palubeckis2004bMST5(const QUBOInstance& qi, double runtime_limit,
		      bool validation, QUBOCallback *qc);
};

class Palubeckis2004bSTS : public QUBOHeuristic {
 public:
  Palubeckis2004bSTS(const QUBOInstance& qi, double runtime_limit,
		     bool validation, QUBOCallback *qc);
};

#endif
