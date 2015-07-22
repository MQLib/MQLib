#include <iostream>
#include "problem/instance.h"
#include "problem/max_cut_instance.h"
#include "problem/qubo_instance.h"

// Load input matrix from provided file
QUBOInstance::QUBOInstance(const std::string& filename) {
  Instance::Load(filename, &nonzero_, &all_nonzero_, &lin_, false);
}

// Load input matrix from provided on- and off-diagonal entries
QUBOInstance::QUBOInstance(const std::vector<Instance::InstanceTuple>& offDiag,
                           const std::vector<double>& mainDiag, int dimension) :
  lin_(mainDiag) {
  if (mainDiag.size() != dimension) {
    std::cout << "Illegal dimension on main diagonal vector" << std::endl;
    exit(1);
  }
  Instance::Load(dimension, offDiag, &nonzero_, &all_nonzero_, NULL, false);
}


QUBOInstance::QUBOInstance(const MaxCutInstance& mi) {
  // Initialize class data structures based on size of base graph
  for (int count=0; count < mi.get_size(); ++count) {
    nonzero_.push_back(std::vector<std::pair<int, double> >());
    lin_.push_back(0.0);
  }

  // Process all edges in input graph
  for (auto iter = mi.get_all_edges_begin(); iter != mi.get_all_edges_end(); 
       ++iter) {
    AddMaxCutEdge(iter->first.first, iter->first.second, iter->second);
  }
}

void QUBOInstance::AddMaxCutEdge(int i, int j, double w_ij) {
  // Every maxcut edge adds its weight to the linear term of the corresponding
  // endpoint variables, but has negative its weight as the quadratic
  // interaction term between the variables. For details on this construction,
  // see equivalence.pdf.
  if (i != j) {
    lin_[i] += w_ij;
    lin_[j] += w_ij;
    nonzero_[i].push_back(std::pair<int, double>(j, -1.0 * w_ij));
    nonzero_[j].push_back(std::pair<int, double>(i, -1.0 * w_ij));
    all_nonzero_.push_back(std::pair<std::pair<int, int>, double>(std::pair<int, int>(i, j), -1.0 * w_ij));
  }
}

// Copy constructor
QUBOInstance::QUBOInstance(const QUBOInstance& qi) :
  nonzero_(qi.nonzero_),
  all_nonzero_(qi.all_nonzero_),
  lin_(qi.lin_) {}

// Copy assignment constructor
QUBOInstance& QUBOInstance::operator=(const QUBOInstance& qi) {
  nonzero_ = qi.nonzero_;
  all_nonzero_ = qi.all_nonzero_;
  lin_ = qi.lin_;
  return *this;
}
