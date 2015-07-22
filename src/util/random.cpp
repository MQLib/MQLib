#include <math.h>
#include <algorithm>
#include <vector>
#include "util/random.h"

int Random::RouletteWheel(const std::vector<double>& scores) {
  double score_sum = 0.0;
  for (auto it = scores.begin(); it != scores.end(); it++) {
    score_sum += *it;
  }
  double r = Random::RandDouble();
  double cum_sum = 0.0;
  for (size_t i = 0; i < scores.size(); i++) {
    cum_sum += scores[i] / score_sum;
    if (r <= cum_sum) {
      return i;
    }
  }
  return scores.size() - 1;
}

void Random::MultiRouletteWheel(const std::vector<double>& scores, int m,
				std::vector<int>* indices) {
  // From comment by Dinre at http://stackoverflow.com/questions/15113650/faster-weighted-sampling-without-replacement; the Efraimidis & Spirakis paper.
  // Paper: http://www.sciencedirect.com/science/article/pii/S002001900500298X#

  // The first value is U[0, 1]^(1/weight), which is the score assigned to
  // a particular element. The second element is the list index. We will not
  // compute scores for non-positive elements (we'll just remove them), making
  // the need for the index. We want to select the largest m elements; for
  // convenience, we'll actually negate the first term and use the default
  // behavior, which is the smallest m elements.

  // For numerical stability, we'll actually use ln(U[0, 1]^(1/weight)) (which
  // is monotone in the original score). That is equal to ln(U[0, 1]) / weight.
  // Again, we'll negate this.
  std::vector<std::pair<double, int> > values;
  for (int ct=0; ct < scores.size(); ++ct) {
    if (scores[ct] > 0.0) {
      values.push_back(std::pair<double, int>(-log(RandDouble()) / scores[ct],
					      ct));
    }
  }
  
  // partial_sort will make the first m elements be the smallest. Only run if
  // m is smaller than the total number of elements. Runs in n lg(m) time, which
  // is better than n lg(n) for just sorting.
  if (m < values.size()) {
    std::partial_sort(values.begin(), values.begin() + m, values.end());
  }

  // Return the first m elements (or all of them, if there are not enough)
  int num_return = std::min<int>(m, values.size());
  indices->clear();
  for (int ct=0; ct < num_return; ++ct) {
    indices->push_back(values[ct].second);
  }
}
