#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>

#include "util/randomForest.h"

RandomForest::RandomForest(const std::string& filename) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    std::cout << "File cannot be opened: " << filename << std::endl;
    exit(1);    
  }

  std::string line;
  int nnodes;

  // First line
  // Element 0: number of random forest trees
  // Element 1: number of variables
  // Element 2: number of nodes across all random forest trees
  getline(file, line);
  if (sscanf(line.c_str(), "%d %d %d", &_ntree, &_nvar, &nnodes) != 3) {
    std::cout << "Illegal first line in random forest file: " << line <<
      std::endl;
    exit(1);
  }

  // Second line: space-separated offsets of each tree's start within the nodes
  getline(file, line);
  std::istringstream is(line);
  _offset = std::vector<int>(std::istream_iterator<int>(is),
                             std::istream_iterator<int>());
  if (_offset.size() != _ntree) {
    std::cout << "Wrong number of tree offsets in " << filename << std::endl;
    exit(1);
  }

  // Subsequent lines (one for each node)
  // Element 0: left daughter (0-indexed; -1 for none)
  // Element 1: right daughter (0-indexed; -1 for none)
  // Element 2: split variable (0-indexed; -1 for predict 0; -2 for predict 1)
  // Element 3: split point
  for (int ct=0; ct < nnodes; ++ct) {
    int left, right, var;
    double split;
    getline(file, line);
    if (sscanf(line.c_str(), "%d %d %d %lf", &left, &right, &var, &split) != 4) {
      std::cout << "Illegal node line: " << line << std::endl;
      exit(1);
    }
    _left.push_back((short)left);
    _right.push_back((short)right);
    _var.push_back((short)var);
    _split.push_back(split);
  }
}

double RandomForest::Predict(const std::vector<double>& vars) {
  if (vars.size() != _nvar) {
    std::cout << "Wrong number of variables in RandomForest::Predict" <<
      std::endl;
    exit(1);
  }

  // Number of trees predicting positive;
  int positive = 0;

  // Determine prediction for each tree
  for (int tree=0; tree < _ntree; ++tree) {
    int offset = _offset[tree];
    int line = offset;
    int var = _var[line];
    while (var >= 0) {
      if (vars[var] <= _split[line]) {
        // Move to left neighbor
        line = offset + _left[line];
      } else {
        // Move to right neighbor
        line = offset + _right[line];
      }
      var = _var[line];
    }

    // var contains -1 for negative prediction and -2 for positive prediction
    if (var == -2) {
      ++positive;
    }
  }

  // Final prediction is proportion of trees predicting positive
  return ((double)positive) / _ntree;
}
