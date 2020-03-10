#include <string.h>
#include <iostream>
#include <vector>
#include "heuristics/maxcut/duarte2005.h"
#include "util/random.h"

///////////////////////////////////////////////////////////////////////////////
// SOLUTION
Duarte2005Solution::Duarte2005Solution(const MaxCutSolution &x) :
  MaxCutSolution(x) {}

// PAPER: procedure Local_Search(g). Note that this differs from
// ExtendedSolution::AllFirst1Swap() by only making one pass through the nodes.
void Duarte2005Solution::Greedy1Swap() {
  // Note:  1 -> \in C,  -1 -> \in C'
  for (int v = 0; v < N_; v++) {
    // We already have diff for every swap pre-calculated
    // PAPER: sigma (v) = \sum_(u \in C ) w_{v,u}
    //        sigma'(v) = \sum_(u \in C') w_{v,u}
    // Is there benefit of moving between sets?
    if (ImprovingMove(v)) {
      UpdateCutValues(v);
    }
  }
}

// PAPER: Child=FixCross(Father,Mother)
Duarte2005Solution::Duarte2005Solution(const Duarte2005Solution& father,
                                       const Duarte2005Solution& mother) :
  MaxCutSolution(father) {
  for (int v = 0; v < N_; ++v) {
    // PAPER: f(A,B) = { 0  if A = 0 and B = 0
    //                 { 1  if A = 1 and B = 1
    //                 { Coin toss otherwise
    if (mother.assignments_[v] != father.assignments_[v] &&
	Random::RandDouble() < 0.5) {
      // Mother and father differ; according to coin flip, reverse father value
      UpdateCutValues(v);
    }
  }
}

// Mutation operator
void Duarte2005Solution::Mutate(double flip_chance) {
  for (int v = 0; v < mi_.get_size(); v++)
    if (Random::RandDouble() < flip_chance)
      UpdateCutValues(v);
}


// VNS
void Duarte2005Solution::VNS(int k_max) {
  // PAPER: /*First Neighbourhood Strucutre*/
  //        k = 1;
  int k = 1;
  // PAPER: while k < kmax do
  while (k < k_max) {
    // PAPER: /*Select an random solution in k-
    //        neighbourhood structure*/
    //        x' = Random(x,Nk(x))
    Duarte2005Solution x_dash(*this);
    // PAPER: In the case of the Max-Cut problem, the kth-order
    //        neighbourhood is defined by all solutions that can be
    //        derived from the current one by selecting k vertices
    //        and transferring each vertex from one subset of the
    //        vertex bipartition to the other subset
    // Succinctly: pick k indices randomly and flip them
    std::vector<bool> flipped(mi_.get_size(), false);
    for (int flip = 0; flip < k; flip++) {
      int idx = Random::RandInt(0, mi_.get_size()-1);
      if (!flipped[idx]) {
        x_dash.UpdateCutValues(idx);
      }
    }
    // PAPER: Use the local search procedure shown in Figure 3
    //        x'' = LocalSearch(x')
    x_dash.Greedy1Swap();

    // Take it if its better
    if (x_dash.ImprovesOver(*this)) {
      operator=(x_dash);  // Copy x_dash over the current solution
      k = 1;
    } else {
      k = k + 1;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// HEURISTIC

// Duarte2005 is an "evolutionary algorithm based on low-level hybridization
// between a memetic algorithm and Variable Neighbourhood Search"
Duarte2005::Duarte2005(const MaxCutInstance& mi, double runtime_limit,
		       bool validation, MaxCutCallback *mc) :
  MaxCutHeuristic(mi, runtime_limit, validation, mc) {
  // Parameters as given in Section 8
  // Initial population size
  const int init_pop_size = 50;
  // Local Search Probability
  const double p_i = 0.25;
  // Probability of crossover
  const double p_crossover = 0.6;
  // Maximum number of generations
  const int max_generations = 50;
  // Probability of mutation process
  const double p_m = 1.0 / mi.get_size();
  // Maximum neighbourhood order
  const int k_max = 1 + mi.get_size() / 100;

  // Repeat until termination criterion is met
  for (int iter=0; MaxCutHeuristic::Report(iter); ++iter) {
    // Generate random initial population
    // PAPER: gg=Intial_Population();
    std::vector<Duarte2005Solution> population;
    for (int sol = 0; sol < init_pop_size; sol++) {
      population.push_back(Duarte2005Solution::RandomSolution(mi, this));
    }
    
    // Improve a pi fraction of the initial population members with local search
    // PAPER: /*Optimize initial population*/
    //        Apply(Local_Search(), pi)
    //        Evaluate_Population();
    //        Best_Solution = Best_Individual();
    std::vector<double> population_scores(init_pop_size, 0.0);
    double best_score = 0;
    for (int sol = 0; sol < init_pop_size; sol++) {
      if (Random::RandDouble() < p_i) {
	population[sol].Greedy1Swap();
      }
      population_scores[sol] = population[sol].get_weight();
      if (population_scores[sol] > best_score) {
	best_score = population_scores[sol];
      }
    }

    // PAPER: for i = 1 to MaxGen
    for (int gen = 0; gen < max_generations; gen++) {
      // With probability p_crossover add a crossed-over solution to the new
      // population and with probability 1-p_crossover keep a randomly
      // selected element
      std::vector<Duarte2005Solution> next_population;
      for (int i=0; i < init_pop_size; ++i) {
	// PAPER: /*Criteria: Random Wheel*/
	int father = Random::RouletteWheel(population_scores);
	// PAPER: r = rand01();/*Random function*/
	double r = Random::RandDouble();
	// PAPER: if (r < pc)
	if (r < p_crossover) {
	  // PAPER: Criteria: Random Wheel
	  int mother = Random::RouletteWheel(population_scores);
	  // PAPER: Child=FixCross(Father,Mother)
	  Duarte2005Solution child =
	    Duarte2005Solution::FixCross(population[father],
					 population[mother]);
	  // PAPER: Apply(VNS(Child, pi))
	  if (Random::RandDouble() < p_i) {
	    child.VNS(k_max);
	  }
	  // PAPER: InsertInPopulation(Child)
	  next_population.push_back(child);
	} else {
	  // PAPER: InsertInPopulation(Father)
	  next_population.push_back(population[father]);
	}  // PAPER: end if
      }  // PAPER: end while
      
      // Though the paper mutates before updating the best solution, this has
      // a good chance of mutating away from a good solution, so we'll
      // evaluate the population here.
      // PAPER: Evaluate_Population()
      //        Best_Solution = Best_Invidual
      population = next_population;
      population_scores = std::vector<double>(population.size());
      int best_score_sol = 0;
      best_score = 0;
      for (int sol = 0; sol < population.size(); sol++) {
	const double score_sol = population[sol].get_weight();
	population_scores[sol] = score_sol;
	if (score_sol > best_score) {
	  best_score = score_sol;
	  best_score_sol = sol;
	}
      }
      
      // Check if termination criterion met
      if (!Report(population[best_score_sol], iter)) {
        return;
      }
      
      // PAPER: Apply(Mutation(),pm)
      for (int sol = 0; sol < next_population.size(); sol++)
	next_population[sol].Mutate(p_m);
    } // PAPER: end for
  }
}
