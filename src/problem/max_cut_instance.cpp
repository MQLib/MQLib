#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include "problem/instance.h"
#include "problem/max_cut_instance.h"
#include "problem/qubo_instance.h"

// Load instance from file
MaxCutInstance::MaxCutInstance(const std::string& filename) {
  Instance::Load(filename, &edges_, &all_edges_, NULL, false);
}

// Load instance from edge list and dimension
MaxCutInstance::MaxCutInstance(const std::vector<Instance::InstanceTuple>& edgeList,
                               int dimension) {
  Instance::Load(dimension, edgeList, &edges_, &all_edges_, NULL, false);
}

// Convert QUBOInstance to MaxCutInstance
MaxCutInstance::MaxCutInstance(const QUBOInstance& qi) {
  // Initialize class data structures based on size of qi and linear values
  std::vector<double> masterNodeWeights(qi.get_lin());  // Weights to added node
  for (int count=0; count < qi.get_size() + 1; ++count) {
    edges_.push_back(std::vector<std::pair<int, double> >());
  }

  // Process all non-zero entries of qi
  for (auto iter = qi.get_all_nonzero_begin(); iter != qi.get_all_nonzero_end();
       ++iter) {
    AddQUBONonzero(iter->first.first, iter->first.second, iter->second,
		   &masterNodeWeights);
  }

  // Add edges to master node (the last node)
  int master = qi.get_size();
  for (int count=0; count < qi.get_size(); ++count) {
    if (masterNodeWeights[count] != 0.0) {
      edges_[count].push_back(std::pair<int, double>(master,
						     masterNodeWeights[count]));
      edges_[master].push_back(std::pair<int, double>(count,
						      masterNodeWeights[count]));
      all_edges_.push_back(std::pair<std::pair<int, int>, double>(std::pair<int, int>(count, master), masterNodeWeights[count]));
    }
  }
}

void MaxCutInstance::AddQUBONonzero(int i, int j, double q_ij,
				    std::vector<double>* masterNodeWeights) {
  (*masterNodeWeights)[i] += q_ij;
  (*masterNodeWeights)[j] += q_ij;
  edges_[i].push_back(std::pair<int, double>(j, -1.0 * q_ij));
  edges_[j].push_back(std::pair<int, double>(i, -1.0 * q_ij));
  all_edges_.push_back(std::pair<std::pair<int, int>, double>(std::pair<int, int>(i, j), -1.0 * q_ij));
}

// Shuffling the edge sets
void MaxCutInstance::GetShuffledEdges(std::vector<std::pair<std::pair<int, int>, double> >* ret) const {
  *ret = all_edges_;
  random_shuffle(ret->begin(), ret->end());
}

bool SortCompare(const std::pair<std::pair<int, int>, double>& i,
		 const std::pair<std::pair<int, int>, double>& j) {
  return (i.second > j.second);
}

void MaxCutInstance::GetSortedEdges(std::vector<std::pair<std::pair<int, int>, double> >* ret) const {
  *ret = all_edges_;
  sort(ret->begin(), ret->end(), SortCompare);
}

// Copy constructor
MaxCutInstance::MaxCutInstance(const MaxCutInstance& mi) :
  edges_(mi.edges_),
  all_edges_(mi.all_edges_) {}

// Copy assignment constructor
MaxCutInstance& MaxCutInstance::operator=(const MaxCutInstance& mi) {
  edges_ = mi.edges_;
  all_edges_ = mi.all_edges_;
  return *this;
}

bool MaxCutInstance::CheckGraph() const {
  int n = get_size();
  std::unordered_set<int> edges;  // Hash edge (i, j), i < j, as n*i+j
  for (auto iter=all_edges_.begin(); iter != all_edges_.end(); ++iter) {
    int i = std::min<int>(iter->first.first, iter->first.second);
    int j = std::max<int>(iter->first.first, iter->first.second);
    std::pair<std::unordered_set<int>::iterator, bool> result =
      edges.insert(i*n + j);
    if (!result.second) {
      // Element was already in set
      std::cout << "Repeated edge: " << iter->first.first << " -> " <<
        iter->first.second << std::endl;
      return false;
    }
  }
  return true;
}

void MaxCutInstance::PrintInstance() const {
  // Number of nodes and edges
  std::cout << edges_.size() << " " << all_edges_.size() << std::endl;
  std::cout.precision(15);
  for (auto iter=all_edges_.begin(); iter != all_edges_.end(); ++iter) {
    std::cout << iter->first.first+1 << " " << iter->first.second+1 << " " <<
      iter->second << std::endl;
  }
}
