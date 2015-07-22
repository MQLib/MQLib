#ifndef UTIL_RANDOM_FOREST_H_
#define UTIL_RANDOM_FOREST_H_

#include <string>
#include <vector>

class RandomForest {
 public:
  // Load a random forest from a file; the random forest should have been saved
  // with the store.rf R function in the scripts folder.
  RandomForest(const std::string& filename);

  // Predict the probability of true from the random forest given the passed
  // independent variable values.
  double Predict(const std::vector<double>& vars);

 private:
  // Number of trees in random forest
  int _ntree;

  // Number of variables in prediction problem
  int _nvar;

  // Offset of each tree's start in the _left, _right, _var, and _split vectors
  std::vector<int> _offset;

  // Left neighbor of each node
  std::vector<short> _left;

  // Right neighbor of each node
  std::vector<short> _right;

  // Split variable of each node (-1: predict negative; -2: predict positive)
  std::vector<short> _var;

  // Split value of each node
  std::vector<double> _split;
};

#endif
