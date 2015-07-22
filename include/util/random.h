#ifndef UTIL_RANDOM_H_
#define UTIL_RANDOM_H_

#include <stdlib.h>

class Random {
 public:
  // Return a random double in [0, 1)
  static inline double RandDouble() {
    return ((double)rand()) / (((long)RAND_MAX)+1);
  }
  
  // Return a random double in [min_val, max_val)
  static inline double RandDouble(double min_val, double max_val) {
    return min_val + (max_val-min_val) * RandDouble();
  }
  
  // Return a random int in [min_val, max_val]
  static inline int RandInt(int min_val, int max_val) {
    return min_val + (rand() % (max_val - min_val + 1));
  }

  // Roulette Wheel Selection -- the passed vector is the weight of each
  // selection; this makes a weighted random selection and returns the index.
  static int RouletteWheel(const std::vector<double>& scores);
  
  // Multi roulette wheel selection -- select m elements without replacement,
  // with each selection weighted by the scores vector. Return the selected
  // indices in the passed indices vector.
  static void MultiRouletteWheel(const std::vector<double>& scores, int m,
				 std::vector<int>* indices);

};

#endif
