#include <math.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include "heuristics/qubo/palubeckis2006.h"
#include "util/random.h"

Palubeckis2006Solution::Palubeckis2006Solution(const QUBOSolution &x) :
  Palubeckis2004bSolution(x) {}

void Palubeckis2006Solution::GSP(int rTilde, int b) {
  std::vector<bool> I(N_, true);  // Whether each var is included in set I
  for (int r=0; r < rTilde; ++r) {
    // Step 2: Sort elements of I by their diff_weight_ values (we will be
    // selecting the b largest values). Because we only want the largest b
    // values, we will use partial_sort, which should run in n lg(b) time
    // instead of the n lg(n) we would get from just sorting. Use effective_b
    // to limit to the number of nodes we have, and exit if no nodes have I set.
    std::vector<std::pair<double, int> > sorted;
    for (int i=0; i < N_; ++i) {
      if (I[i]) {
	sorted.push_back(std::pair<double, int>(diff_weights_[i], i));
      }
    }
    if (sorted.size() == 0) {
      break;
    }
    int effective_b = std::min<int>(b, sorted.size());
    std::partial_sort(sorted.begin(), sorted.begin()+effective_b, sorted.end(),
		      std::greater<std::pair<double, int> >());

    // Step 3: Randomly select q from J (the first b elements of sorted). Flip
    // its set inclusion and remove q from I.
    int q = sorted[Random::RandInt(0, effective_b-1)].second;
    UpdateCutValues(q);
    I[q] = false;
  }
}

// ITS method from Palubeckis2006
Palubeckis2006::Palubeckis2006(const QUBOInstance& qi, double runtime_limit,
			       bool validation, QUBOCallback *qc) :
  QUBOHeuristic(qi, runtime_limit, validation, qc) {
  // Step 0: Setup parameters
  int mu;
  if (qi.get_size() <= 3000) {
    mu = 10000;
  } else if (qi.get_size() <= 5000) {
    mu = 12000;
  } else {
    mu = 15000;
  }
  int mtilde = mu * qi.get_size();
  int d1 = 10;
  double d2 = 0.1;
  int rangemax = std::max((int)floor(d2*qi.get_size()), d1);
  int b = 5;

  // Since ITS is perturb-then-optimize (GSP to perturb and TS to optimize), we
  // will wrap it in a random restart. The termination criterion used in the
  // paper is based on runtime, making it tricky to exactly replicate. We'll do
  // this in two steps: 1) We'll compute a runtime limit for the hardware in the
  // paper based on the problem instance size, and 2) We'll convert this limit
  // to our hardware using an (approximate) conversion from Dongarra (2014).
  
  // Part 1: Compute a runtime limit for Palubeckis2006 hardware. The runtime
  // limits reported in the results section were 600, 900, 1800, 3600, 5400, and
  // 9000 seconds, respectively, for instance node counts 2500, 3000, 4000,
  // 5000, 6000, and 7000. This is mildly quadratic, and we can use linear
  // regression with n and n^2 terms to get runtime limit equation:
  // limit = 2490.8 - 1.635305n + 0.0003633706n^2. This model does a fine job
  // of interpolating the runtime limits (adj. R^2 = 0.9935) and seems reasonable
  // for n larger than 7000, but it has a minimum at 2250 and increases for
  // smaller values of n, which doesn't make any sense. To address this, we'll
  // use linear interpolation of the runtime limit for n in the range [0, 2500],
  // with limit 1 second for the n=0 case. This yields equation:
  // limit = 1 + 0.2690411n.
  double limit_raw;
  if (qi.get_size() < 2500) {
    limit_raw = 1.0 + 0.2690411 * qi.get_size();
  } else {
    limit_raw = 2490.8 - 1.635305 * qi.get_size() +
      0.0003633706 * qi.get_size() * qi.get_size();
  }

  // Part 2: Covert this limit to our hardware using an (approximate) conversion
  // from Dongarra (2014). The reported hardware is "a Pentium III 800 PC"; we
  // don't have access to such a machine so we'll instead use reported flop
  // counts from Dongarra (2014)'s LINPACK benchmarks. While the benchmarks
  // don't include this exact model, they include similar computers (flop counts
  // are for the LINPACK 100 benchmark):
  // Intel Pentium III 933 MHz (234 Mflop/s)
  // Intel Pentium III 933 MHz (192 Mflop/s)
  // Pentium III (750 MHz) (138 Mflop/s)
  // Dell Dimension XPS T500 (500 MHz) Pentium III (86 Mflop/s)
  // Intel Pentium III 550 MHz (86 Mflop/s)
  // Intel Pentium III 550 MHz (80 Mflop/s)
  // Intel Pentium III 550 MHz (74 Mflop/s)
  // Dell Pentium III 450 MHz (65 Mflop/s)
  // Using linear interpolation on the clock speeds (adj. R^2 = 0.9318), we
  // estimate the LINPACK 100 benchmark value for a Pentium III 800 PC to be
  // 165.8905 Mflop/s. The m3.medium nodes used in our computational experiments
  // have LINPACK 100 performance 1128.9522 Mflop/s. We can use this information
  // to obtain an adjusted runtime limit for use in random restarts:
  double limit = limit_raw * 165.8905 / 1128.9522;

  while (true) {
    double start = Runtime();

    // Step 1: Randomly generate solution
    Palubeckis2006Solution x = QUBOSolution::RandomSolution(qi, this);
    double best_objective = x.get_weight();

    while (Runtime() - start <= limit) {
      // Step 3: Call TS on solution
      x.TS(&best_objective, mtilde);
      
      // Step 4: Check if stepping criterion is satisfied
      if (!Report(x)) {
        return;
      }
      
      // Step 5: Call GSP
      x.GSP(Random::RandInt(d1, rangemax), b);
    }
  }
}
