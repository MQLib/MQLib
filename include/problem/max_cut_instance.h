#ifndef PROBLEM_MAX_CUT_INSTANCE_H_
#define PROBLEM_MAX_CUT_INSTANCE_H_

#include <string>
#include <utility>
#include <vector>
#include "problem/instance.h"

// Forward declaration (since MaxCutInstance can be constructed from QUBOInstance
// and vice versa)
class QUBOInstance;

class MaxCutInstance {
 public:
  // Constructor (load graph from passed file)
  MaxCutInstance(const std::string& filename);

  // Constructor (provide graph by providing edge list through vector of tuples)
  MaxCutInstance(const std::vector<Instance::InstanceTuple>& edgeList,
                 int dimension);

  // Construct from a QUBO instance
  MaxCutInstance(const QUBOInstance& qi);

  // Copy constructor
  MaxCutInstance(const MaxCutInstance &mi);

  // Copy assignment constructor
  MaxCutInstance& operator=(const MaxCutInstance& mi);

  // Populate the passed vector with all the edges, and shuffle.
  void GetShuffledEdges(std::vector<std::pair<std::pair<int, int>, double> >*
			to_shuffle) const;
  void GetSortedEdges(std::vector<std::pair<std::pair<int, int>, double> >*
		      to_sort) const;

  // Check if the graph has any repeated edges
  bool CheckGraph() const;

  // Output the adjacency list (with or without weights) of a graph to the
  // screen
  void PrintInstance() const;

  // Getters
  int get_size() const {  return edges_.size();  }
  int get_edge_count() const {  return all_edges_.size();  }
  int get_vertex_degree(int idx) const {  return edges_[idx].size();  }

  std::vector<std::vector<std::pair<int, double> > > get_edges() const { return edges_; }
  std::vector<std::pair<std::pair<int, int>, double> > get_all_edges() const {return all_edges_; }


  std::vector<std::pair<int, double> >::const_iterator get_edges_begin(int idx)
    const {
    return edges_[idx].begin();
  }
  std::vector<std::pair<int, double> >::const_iterator get_edges_end(int idx)
    const {
    return edges_[idx].end();
  }  

  std::vector<std::pair<std::pair<int, int>, double> >::const_iterator
    get_all_edges_begin() const {  return all_edges_.begin();  }

  std::vector<std::pair<std::pair<int, int>, double> >::const_iterator
    get_all_edges_end() const {  return all_edges_.end();  }

 protected:
  // During construction from a QUBO instance, add non-zero matrix value q_ij,
  // updating the edge weights as well as the weight to the added "master node."
  void AddQUBONonzero(int i, int j, double q_ij,
		      std::vector<double>* masterNodeWeights);

  // Non-const iterators
  std::vector<std::pair<int, double> >::iterator
    get_mutable_edges_begin(int idx) {  return edges_[idx].begin();  }
  std::vector<std::pair<int, double> >::iterator
    get_mutable_edges_end(int idx) {  return edges_[idx].end();  }
  std::vector<std::pair<std::pair<int, int>, double> >::iterator
    get_mutable_all_edges_begin() {  return all_edges_.begin();  }

  std::vector<std::pair<std::pair<int, int>, double> >::iterator
    get_mutable_all_edges_end() {  return all_edges_.end();  }

  std::vector<std::vector<std::pair<int, double> > > edges_;
  std::vector<std::pair<std::pair<int, int>, double> > all_edges_;
};

#endif
