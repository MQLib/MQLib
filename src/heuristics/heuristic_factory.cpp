#include <algorithm>
#include <iostream>
#include <vector>

#include "heuristics/heuristic_factory.h"

#include "heuristics/maxcut/baseline.h"
#include "heuristics/maxcut/burer2002.h"
#include "heuristics/maxcut/deSousa2013.h"
#include "heuristics/maxcut/duarte2005.h"
#include "heuristics/maxcut/festa2002.h"
#include "heuristics/maxcut/laguna2009.h"

#include "heuristics/qubo/alkhamis1998.h"
#include "heuristics/qubo/beasley1998.h"
#include "heuristics/qubo/glover1998a.h"
#include "heuristics/qubo/glover2010.h"
#include "heuristics/qubo/hasan2000.h"
#include "heuristics/qubo/katayama2000.h"
#include "heuristics/qubo/katayama2001.h"
#include "heuristics/qubo/lodi1999.h"
#include "heuristics/qubo/lu2010.h"
#include "heuristics/qubo/merz1999.h"
#include "heuristics/qubo/merz2002.h"
#include "heuristics/qubo/merz2004.h"
#include "heuristics/qubo/palubeckis2004b.h"
#include "heuristics/qubo/palubeckis2006.h"
#include "heuristics/qubo/pardalos2008.h"

// Get a pointer to a heuristic of a templated type
template<typename T> MaxCutHeuristic* NewMaxCutHeur(const MaxCutInstance& mi,
                                                    double runtime_limit,
                                                    bool validation,
                                                    MaxCutCallback *mc)  { 
  return new T(mi, runtime_limit, validation, mc);
}
template<typename T> QUBOHeuristic* NewQUBOHeur(const QUBOInstance& qi,
                                                double runtime_limit,
                                                bool validation,
                                                QUBOCallback *qc)  {
  return new T(qi, runtime_limit, validation, qc);
}

HeuristicFactory::HeuristicFactory() {
  // Register all Max-Cut heuristics
  max_cut_map_["BURER2002"] =
    MaxCutCreator(&NewMaxCutHeur<Burer2002>,
                  "Rank-2 relaxation of G-W algorithm, with local search");
  max_cut_map_["FESTA2002G"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002G>, "GRASP with local search");
  max_cut_map_["FESTA2002GPR"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002GPR>, "GRASP with path-relinking");
  max_cut_map_["FESTA2002VNS"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002VNS>, "VNS");
  max_cut_map_["FESTA2002VNSPR"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002VNSPR>, "VNS with path-relinking");
  max_cut_map_["FESTA2002GVNS"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002GVNS>, "GRASP with VNS local search");
  max_cut_map_["FESTA2002GVNSPR"] =
    MaxCutCreator(&NewMaxCutHeur<Festa2002GVNSPR>,
                  "GRASP & VNS with path-relinking");
  max_cut_map_["DUARTE2005"] =
    MaxCutCreator(&NewMaxCutHeur<Duarte2005>,
                  "Genetic algorithm with VNS as local search");
  max_cut_map_["LAGUNA2009CE"] =
    MaxCutCreator(&NewMaxCutHeur<Laguna2009CE>, "Cross-entropy method");
  max_cut_map_["LAGUNA2009HCE"] =
    MaxCutCreator(&NewMaxCutHeur<Laguna2009HCE>,
                  "Cross-entropy method with local search");
  max_cut_map_["DESOUSA2013"] =
    MaxCutCreator(&NewMaxCutHeur<deSousa2013>,
                  "Estimation of distribution algorithm");
  max_cut_map_["BASELINE"] =
    MaxCutCreator(&NewMaxCutHeur<Baseline>, "Baseline heuristic");

  // Register all QUBO heuristics
  qubo_map_["PALUBECKIS2004bMST1"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bMST1>, "Tabu search procedure");
  qubo_map_["PALUBECKIS2004bMST2"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bMST2>, "Iterated tabu search");
  qubo_map_["PALUBECKIS2004bMST3"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bMST3>, "Tabu search with GRASP");
  qubo_map_["PALUBECKIS2004bMST4"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bMST4>,
                "Tabu search with long-term memory");
  qubo_map_["PALUBECKIS2004bMST5"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bMST5>, "Iterated tabu search");
  qubo_map_["PALUBECKIS2004bSTS"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2004bSTS>, "Tabu search procedure");
  qubo_map_["MERZ1999GLS"] =
    QUBOCreator(&NewQUBOHeur<Merz1999GLS>,
                "Genetic algorithm, with crossover and local search");
  qubo_map_["MERZ1999MUTATE"] =
    QUBOCreator(&NewQUBOHeur<Merz1999Mutation>,
                "Genetic algorithm, with mutation only");
  qubo_map_["MERZ1999CROSS"] =
    QUBOCreator(&NewQUBOHeur<Merz1999Crossover>,
                "Genetic algorithm, with crossover only");
  qubo_map_["PALUBECKIS2006"] =
    QUBOCreator(&NewQUBOHeur<Palubeckis2006>, "Iterated tabu search");
  qubo_map_["KATAYAMA2001"] =
    QUBOCreator(&NewQUBOHeur<Katayama2001>, "Simulated annealing");
  qubo_map_["LU2010"] =
    QUBOCreator(&NewQUBOHeur<Lu2010>, "Genetic algorithm with tabu search");
  qubo_map_["PARDALOS2008"] =
    QUBOCreator(&NewQUBOHeur<Pardalos2008>, "Global equilibrium search");
  qubo_map_["GLOVER2010"] =
    QUBOCreator(&NewQUBOHeur<Glover2010>, "Tabu search with long-term memory");
  qubo_map_["BEASLEY1998SA"] =
    QUBOCreator(&NewQUBOHeur<Beasley1998SA>, "Simulated annealing");
  qubo_map_["BEASLEY1998TS"] =
    QUBOCreator(&NewQUBOHeur<Beasley1998TS>, "Tabu search");
  qubo_map_["MERZ2002GREEDY"] =
    QUBOCreator(&NewQUBOHeur<Merz2002Greedy>, "GRASP without local search");
  qubo_map_["MERZ2002ONEOPT"] =
    QUBOCreator(&NewQUBOHeur<Merz2002OneOpt>,
                "1-opt local search with random restarts");
  qubo_map_["MERZ2002KOPT"] =
    QUBOCreator(&NewQUBOHeur<Merz2002KOpt>,
                "k-opt local search with random restarts");
  qubo_map_["MERZ2002GREEDYKOPT"] =
    QUBOCreator(&NewQUBOHeur<Merz2002GreedyKOpt>,
                "k-opt local search with GRASP");
  qubo_map_["KATAYAMA2000"] =
    QUBOCreator(&NewQUBOHeur<Katayama2000>,
                "Genetic algorithm with k-opt local search");
  qubo_map_["MERZ2004"] =
    QUBOCreator(&NewQUBOHeur<Merz2004>,
                "Genetic algorithm with k-opt local search");
  qubo_map_["ALKHAMIS1998"] =
    QUBOCreator(&NewQUBOHeur<Alkhamis1998>, "Simulated annealing");
  qubo_map_["HASAN2000GA"] =
    QUBOCreator(&NewQUBOHeur<Hasan2000GA>, "Genetic algorithm");
  qubo_map_["HASAN2000TS"] =
    QUBOCreator(&NewQUBOHeur<Hasan2000TS>, "Tabu search");
  qubo_map_["LODI1999"] =
    QUBOCreator(&NewQUBOHeur<Lodi1999>, "Genetic algorithm");
  qubo_map_["GLOVER1998a"] =
    QUBOCreator(&NewQUBOHeur<Glover1998a>, "Tabu search");
}

MaxCutHeuristic* HeuristicFactory::RunMaxCutHeuristic(const std::string& code,
                                                      const MaxCutInstance& mi,
                                                      double runtime_limit,
                                                      bool validation,
                                                      MaxCutCallback *mc) {
  std::map<std::string, MaxCutCreator>::iterator it;
  it = max_cut_map_.find(code);
  if (it != max_cut_map_.end()) {
    return max_cut_map_[code].first(mi, runtime_limit, validation, mc);
  }
  return NULL;  // Heuristic code not found
}

QUBOHeuristic* HeuristicFactory::RunQUBOHeuristic(const std::string& code,
                                                  const QUBOInstance& qi,
                                                  double runtime_limit,
                                                  bool validation,
                                                  QUBOCallback *qc) {
  std::map<std::string, QUBOCreator>::iterator it;
  it = qubo_map_.find(code);
  if (it != qubo_map_.end()) {
    return qubo_map_[code].first(qi, runtime_limit, validation, qc);
  }
  return NULL;  // Heuristic code not found
}

void HeuristicFactory::MaxCutHeuristicCodes(std::vector<std::string>* codes) {
  if (!codes) {
    return;
  }
  codes->clear();
  for (auto iter = max_cut_map_.begin(); iter != max_cut_map_.end(); ++iter) {
    codes->push_back(iter->first);
  }
  std::sort(codes->begin(), codes->end());
}

void HeuristicFactory::QUBOHeuristicCodes(std::vector<std::string>* codes) {
  if (!codes) {
    return;
  }
  codes->clear();
  for (auto iter = qubo_map_.begin(); iter != qubo_map_.end(); ++iter) {
    codes->push_back(iter->first);
  }
  std::sort(codes->begin(), codes->end());
}

void HeuristicFactory::PrintHeuristicCodes() {
  std::cout << "=========================" << std::endl;
  std::cout << "Heuristic codes and descriptions" << std::endl;
  std::cout << "=========================" << std::endl << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Max-Cut heuristics" << std::endl;
  std::cout << "------------------" << std::endl << std::endl;
  std::vector<std::string> max_cut_names;
  for (auto iter = max_cut_map_.begin(); iter != max_cut_map_.end(); ++iter) {
    max_cut_names.push_back(iter->first);
  }
  std::sort(max_cut_names.begin(), max_cut_names.end());
  for (int i=0; i < max_cut_names.size(); ++i) {
    std::cout << max_cut_names[i] << std::endl;
    std::cout << "  " << max_cut_map_[max_cut_names[i]].second << std::endl;
  }
  std::cout << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << "QUBO heuristics" << std::endl;
  std::cout << "---------------" << std::endl << std::endl;
  std::vector<std::string> qubo_names;
  for (auto iter = qubo_map_.begin(); iter != qubo_map_.end(); ++iter) {
    qubo_names.push_back(iter->first);
  }
  std::sort(qubo_names.begin(), qubo_names.end());
  for (int i=0; i < qubo_names.size(); ++i) {
    std::cout << qubo_names[i] << std::endl;
    std::cout << "  " << qubo_map_[qubo_names[i]].second << std::endl;
  }
}
