#include "heuristics/heuristic_factory.h"
#include <string>
#include <vector>
#include <iostream>

using namespace mqlib;

inline void test(const bool& pred)
{
  if (!pred)
    throw std::runtime_error("Test failed.");
}

std::vector<std::string> max_cut_solvers = {
  "BURER2002",
  "FESTA2002G",
  "FESTA2002GPR",
  "FESTA2002VNS",
  "FESTA2002VNSPR",
  "FESTA2002GVNS",
  "FESTA2002GVNSPR",
  "DUARTE2005",
  "LAGUNA2009CE",
  "LAGUNA2009HCE",
  "DESOUSA2013"
};

std::vector<std::string> qubo_solvers = {
  "PALUBECKIS2004bMST1",
  "PALUBECKIS2004bMST2",
  "PALUBECKIS2004bMST3",
  "PALUBECKIS2004bMST4",
  "PALUBECKIS2004bMST5",
  "PALUBECKIS2004bSTS",
  "MERZ1999GLS",
  "MERZ1999MUTATE",
  "MERZ1999CROSS",
  "PALUBECKIS2006",
  "KATAYAMA2001",
  "LU2010",
  "PARDALOS2008",
  "GLOVER2010",
  "BEASLEY1998SA",
  "BEASLEY1998TS",
  "MERZ2002GREEDY",
  "MERZ2002ONEOPT",
  "MERZ2002KOPT",
  "MERZ2002GREEDYKOPT",
  "KATAYAMA2000",
  "MERZ2004",
  "ALKHAMIS1998",
  "HASAN2000GA",
  "HASAN2000TS",
  "LODI1999",
  "GLOVER1998a" 
};

mqlib::HeuristicFactory factory;

void run_max_cut_heuristic_on_max_cut_problem(const std::string& instance, const std::string& heuristic, const double optimal_cost)
{
  std::cout << "solve max cut problem " << instance << " with max cut heuristic " << heuristic << "\n";
  MaxCutHeuristic* heur = factory.RunMaxCutHeuristic(heuristic, MaxCutInstance(instance), 1.0, false, NULL);
  const mqlib::MaxCutSimpleSolution& sol = heur->get_best_solution();
  test(sol.get_weight() <= optimal_cost + 1e-8);
  delete heur;
}

void run_max_cut_heuristic_on_QUBO_problem(const std::string& instance, const std::string& heuristic, const double optimal_cost)
{
  std::cout << "solve QUBO problem " << instance << " with max cut heuristic " << heuristic << "\n";
  MaxCutHeuristic* heur = factory.RunMaxCutHeuristic(heuristic, MaxCutInstance(QUBOInstance(instance)), 1.0, false, NULL);
  const mqlib::MaxCutSimpleSolution& sol = heur->get_best_solution();
  test(sol.get_weight() <= optimal_cost + 1e-8);
  delete heur;
}

void run_QUBO_heuristic_on_max_cut_problem(const std::string& instance, const std::string& heuristic, const double optimal_cost)
{
  std::cout << "solve max cut problem " << instance << " with QUBO heuristic " << heuristic << "\n";
  QUBOHeuristic* heur = factory.RunQUBOHeuristic(heuristic, QUBOInstance(MaxCutInstance(instance)), 1.0, false, NULL);
  const mqlib::QUBOSimpleSolution& sol = heur->get_best_solution();
  test(sol.get_weight() <= optimal_cost + 1e-8);
  delete heur;
}

void run_QUBO_heuristic_on_QUBO_problem(const std::string& instance, const std::string& heuristic, const double optimal_cost)
{
  std::cout << "solve QUBO problem " << instance << " with QUBO heuristic " << heuristic << "\n";
  QUBOHeuristic* heur = factory.RunQUBOHeuristic(heuristic, QUBOInstance(instance), 1.0, false, NULL);
  const mqlib::QUBOSimpleSolution& sol = heur->get_best_solution();
  test(sol.get_weight() <= optimal_cost + 1e-8);
  delete heur;
}

int main()
{
  const std::string max_cut_test_instance = MQLIB_PROJECT_DIR "/bin/sampleMaxCut.txt";
  const double max_cut_test_instance_cost = 3.0;

  const std::string QUBO_test_instance = MQLIB_PROJECT_DIR "/bin/sampleQUBO.txt";
  const double QUBO_test_instance_cost = 6.0;

  for(auto mc_solver : max_cut_solvers) {
    run_max_cut_heuristic_on_max_cut_problem(max_cut_test_instance, mc_solver, max_cut_test_instance_cost);
    run_max_cut_heuristic_on_QUBO_problem(QUBO_test_instance, mc_solver, QUBO_test_instance_cost);
  }

  for(auto qubo_solver : qubo_solvers) {
    run_QUBO_heuristic_on_max_cut_problem(max_cut_test_instance, qubo_solver, max_cut_test_instance_cost);
    run_QUBO_heuristic_on_QUBO_problem(QUBO_test_instance, qubo_solver, QUBO_test_instance_cost);
  } 
}


