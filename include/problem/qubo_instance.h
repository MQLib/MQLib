#ifndef PROBLEM_QUBO_INSTANCE_H_
#define PROBLEM_QUBO_INSTANCE_H_

#include <vector>
#include <string>
#include "problem/instance.h"

// Forward declaration (since MaxCutInstance can be constructed from QUBOInstance
// and vice versa)
class MaxCutInstance;

class QUBOInstance {
 public:
  // Constructor (load input matrix from passed file)
  QUBOInstance(const std::string& filename);

  // Constructor (provide input matrix by providing vector for the main diagonal
  // and tuples for the off-diagonal entries)
  QUBOInstance(const std::vector<Instance::InstanceTuple>& offDiag,
               const std::vector<double>& mainDiag, int dimension);

  // Constructor from a MAX-CUT graph
  QUBOInstance(const MaxCutInstance& mi);

  // Copy constructor
  QUBOInstance(const QUBOInstance &qi);

  // Copy assignment constructor
  QUBOInstance& operator=(const QUBOInstance& qi);

  // Getters
  int get_size() const {  return nonzero_.size();  }
  int get_edge_count() const {  return all_nonzero_.size();  }

  std::vector<std::pair<int, double> >::const_iterator
    get_nonzero_begin(int idx) const {  return nonzero_[idx].begin();  }

  std::vector<std::pair<int, double> >::const_iterator
    get_nonzero_end(int idx) const {  return nonzero_[idx].end();  }

  std::vector<std::pair<std::pair<int, int>, double> >::const_iterator
    get_all_nonzero_begin() const {  return all_nonzero_.begin();  }

  std::vector<std::pair<std::pair<int, int>, double> >::const_iterator
    get_all_nonzero_end() const {  return all_nonzero_.end();  }

  const std::vector<double>& get_lin() const {  return lin_;  }

 protected:
  // During construction from a maxcut instance, add edge from i to j with
  // weight w_ij.
  void AddMaxCutEdge(int i, int j, double w_ij);

  // Non-const iterators
  std::vector<std::pair<int, double> >::iterator
    get_mutable_nonzero_begin(int idx) {  return nonzero_[idx].begin();  }

  std::vector<std::pair<int, double> >::iterator
    get_mutable_nonzero_end(int idx) {  return nonzero_[idx].end();  }
  std::vector<std::pair<std::pair<int, int>, double> >::iterator
    get_mutable_all_nonzero_begin() {  return all_nonzero_.begin();  }

  std::vector<std::pair<std::pair<int, int>, double> >::iterator
    get_mutable_all_nonzero_end() {  return all_nonzero_.end();  }

  std::vector<std::vector<std::pair<int, double> > > nonzero_;
  std::vector<std::pair<std::pair<int, int>, double> > all_nonzero_;
  std::vector<double> lin_;
};

#endif
