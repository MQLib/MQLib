#include <algorithm>
#include <iostream>
#include "heuristics/qubo/merz1999.h"
#include "util/random.h"

Merz1999Solution::Merz1999Solution(const QUBOSolution &x) :
QUBOSolution(x) {}

Merz1999Solution::Merz1999Solution(const QUBOInstance& qi,
                                   const Merz1999Solution& parent_a,
                                   const Merz1999Solution& parent_b,
                                   QUBOHeuristic *heuristic) :
  QUBOSolution(parent_a) {
  // Implements the HUX cross over with restricted local search
  // Store the bits which were identical with parents before commencing local search.
  std::vector<bool> parents_identical(N_, false);
  // http://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)#Uniform_Crossover_and_Half_Uniform_Crossover
  //   In the half uniform crossover scheme (HUX), exactly half of the 
  //   nonmatching bits are swapped. Thus first the Hamming distance (the
  //   number of differing bits) is calculated. This number is divided by two.
  //   The resulting number is how many of the bits that do not match between
  //   the two parents will be swapped.
  int half_hamming_distance = parent_a.SymmetricDifference(parent_b) / 2;
  std::vector<int> genes(N_, 0);
  for (int i = 0; i < N_; i++)
    genes[i] = i;
  // Pick random ordering of genes
  std::random_shuffle(genes.begin(), genes.end());
  int left_to_swap = half_hamming_distance;
  const std::vector<int>& a_genes = parent_a.get_assignments();
  const std::vector<int>& b_genes = parent_b.get_assignments();
  for (int pos = 0; pos < N_; pos++) {
    int i = genes[pos];
    if (a_genes[i] == b_genes[i]) {
      // Parents are the same in this gene
      parents_identical[i] = true;
    } else {
      // Parents are different in this gene
      // The wikipedia page is confusing because it talks about two children
      // But the paper only has one
      // So, by default we'll take the father
      if (left_to_swap > 0) {
	UpdateCutValues(i);
	--left_to_swap;
      }
    }
  }

  // We now have some combination of the parents, now do the local search
  // PAPER:  The local search applied to the resulting offspring after
  //         recombination is restricted to a region of the search space
  //         defined by the two parents: the genes with equal values in 
  //         the two parents are not modified during local search.
  while (1) {
    double best_move = 0.0;
    int best_pos = -1;
    for (int i=0; i < N_; ++i) {
      if (parents_identical[i]) continue;  // Only modification
      if (diff_weights_[i] > best_move) {
        best_move = diff_weights_[i];
        best_pos = i;
      }
    }
    if (best_pos < 0 || !ImprovingMove(best_pos)) {
      // No more profitable moves
      break;
    }
    
    // Update the diff_weights_ variable and objective
    UpdateCutValues(best_pos);
  }
}


Merz1999Solution::Merz1999Solution(const QUBOInstance& qi,
                                   const Merz1999Solution& parent_a,
                                   const Merz1999Solution& parent_b,
                                   QUBOHeuristic *heuristic,
                                   double dummy) :
  QUBOSolution(parent_a) {
  const std::vector<int>& a_genes = parent_a.get_assignments();
  const std::vector<int>& b_genes = parent_b.get_assignments();
  for (int i = 0; i < N_; i++) {
    if (a_genes[i] == b_genes[i]) {
      // Same, do nothing
    } else {
      if (Random::RandDouble() <= 0.5) {
        // Be like the mother, the opposite of what it is currently
        UpdateCutValues(i);
      }
    }
  }
}

void Merz1999Solution::RestartMutate() {
  // Flip a third of the bits at random
  for (int i = 0; i < N_ / 3; i++) {
    int bit_pos = Random::RandInt(0, N_ - 1);
    UpdateCutValues(bit_pos);
  }
}

void Merz1999Solution::Mutate() {
  // Flip a random bit
  int bit_pos = Random::RandInt(0, N_ - 1);
  UpdateCutValues(bit_pos);
}

// version==0 : the local search variant
// version==1 : the crossover variant
// version==2 : the mutation variant
Merz1999::Merz1999(const QUBOInstance& qi, double runtime_limit, bool validation,
		   QUBOCallback *qc, int version) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {

  int POP_SIZE = (version == 0) ? 40 : 100;
  const double RECOMBINATION_RATE = 0.5;

  // PAPER:  begin
  // PAPER:  initialize population P;
  std::vector<Merz1999Solution> population;
  for (int i = 0; i < POP_SIZE; i++)
    population.push_back(QUBOSolution::RandomSolution(qi, this));
  if (version == 0)
    // PAPER:   foreach individual i \in P do i := Local-Search(i);
    for (int i = 0; i < POP_SIZE; i++)
      population[i].AllBest1Swap();  // ExtendedSolution::AllBest1Swap()

  // PAPER:  repeat
  int last_iteration_with_change = 1;
  for (int iteration=0; ; ++iteration) {
    // PAPER:  for i := 1 to #crossovers do
    int num_crossovers = static_cast<int>(population.size() * RECOMBINATION_RATE);
    for (int i = 0; i < num_crossovers; i++) {
      // PAPER:  select two parents i_a, i_b \in P randomly;
      int i_a = Random::RandInt(0, population.size()-1);
      int i_b = i_a;
      while (i_a == i_b)
        i_b = Random::RandInt(0, population.size()-1);
      // Version 0 (local search)
      if (version == 0) {
        // PAPER:  i_c = Crossover(i_a, i_b);
        // PAPER:  Local-Search(i_c);
        // NOTES:  The local search is not the same as the other uses in this paper
        Merz1999Solution cross = Merz1999Solution::HUX(qi, population[i_a], population[i_b], this);
        // PAPER:  add individual i_c to P;
        population.push_back(cross);
      } else if (version == 1) {
        // PAPER:  i_c = Crossover(i_a, i_b);
        Merz1999Solution cross = Merz1999Solution::Crossover(qi, population[i_a], population[i_b], this);
        // PAPER:  add individual i_c to P;
        population.push_back(cross);
      } else if (version == 2) {
        // PAPER:  In our algorithms, we use a simple bit flip
        //         mutation operator that flips a single bit of the bit
        //         string constituting a solution to the BQP.
        // PAPER:  as well as our genetic algorithm solely based on mutation
        // NOTES:  Mutation applied to... one of the parents only, I guess?
        Merz1999Solution child(population[i_a]);
        child.Mutate();
        population.push_back(child);
      }
    // PAPER:  endfor
    }
    // PAPER:  P := select(P);
    // NOTES:  Take best POP_SIZE of current solutions
    // Retrieve (weight, index) pairs
    std::vector<std::pair<double,int>> population_weights(population.size());
    for (int i = 0; i < population.size(); i++) {
      population_weights[i] = std::pair<double,int>(population[i].get_weight(), i);
    }
    // Sort pairs ascending order
    std::sort(population_weights.begin(), population_weights.end());
    // Take last POP_SIZE of list
    std::vector<Merz1999Solution> new_population;
    for (int i = population.size() - 1; i >= population.size()-POP_SIZE; i--) {
      int j = population_weights[i].second;
      new_population.push_back(population[j]);
    }
    if (population[POP_SIZE-1].get_weight() != new_population[POP_SIZE-1].get_weight()) {
      // So the population has changed
      last_iteration_with_change = iteration;
    }
    population = new_population;
    //std::cout << population[0].get_weight() << std::endl;
    //std::cout << population[POP_SIZE-1].get_weight() << std::endl;
   
    // PAPER:  if P converged then
    // NOTES:  this is defined as the average hamming distance of the population
    //         dropping below 10 OR no change in population for more than 30 gens
    double total_hamming_distance = 0.0;
    int total_hamming_calcs = 0;
    for (int i = 0; i < POP_SIZE; i++) {
      for (int j = i+1; j < POP_SIZE; j++) {
        total_hamming_distance += population[i].SymmetricDifference(population[j]);
        total_hamming_calcs++;  // Math is annoying
      }
    }
    //std::cout << "  total_hamming_distance = " << total_hamming_distance << std::endl;
    double hamming_mean = total_hamming_distance / total_hamming_calcs;
    //std::cout << "  hamming_mean = " << hamming_mean << std::endl;
    if (hamming_mean <= 10 || (iteration - last_iteration_with_change) >= 30) {
      // PAPER:  foeach individual i \in P\{best} do i := LocalSearch(Mutate(i))
      for (int i = 1; i < population.size(); i++) {  // Skip best
        // PAPER:  During the restarts, the individuals were mutated 
        //         by flipping a third of all the bits in the bit vector.
        population[i].RestartMutate();
        if (version == 0) population[i].AllBest1Swap();
        if (!Report(population[i])) {
          break;
        }
      }
    }

    // PAPER:  until terminate=true;
    // NOTE: Iterations are very fast with the crossover and mutate versions of the
    //       algorithm, and the best solution never gets worse through time. As a
    //       result, we only report the new best solution every 10 iterations.
    if ((version == 0 || iteration % 10 == 0) && !Report(population[0])) {
      break;
    }
  }
}
