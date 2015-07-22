#include <math.h>
#include <queue>
#include <string.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include "metrics/max_cut_metrics.h"
#include "problem/heuristic.h"
#include "util/random.h"

GraphMetrics::GraphMetrics(const MaxCutInstance& mi) :
  mi_(mi) {}

double GraphMetrics::GetTime(const struct timeval& start) {
  // Runtime in seconds since the passed start time
  struct timeval end;
  gettimeofday(&end, 0);
  return (end.tv_sec - start.tv_sec) + 0.000001 * (end.tv_usec - start.tv_usec);
}

void GraphMetrics::AllMetrics(std::vector<double>* metrics,
                              std::vector<double>* runtimes) {
  srand(0);  // To give consistency between runs of this code
  int n = mi_.get_size();
  double log_n = log(n);
  int m = mi_.get_edge_count();
  double log_m = log(m);
  double avg_degree = ((double)m)/n;

  // Local clustering coefficient statistics
  struct timeval clust_start;
  gettimeofday(&clust_start, 0);
  std::vector<double> clustering;
  GetClusteringData(&clustering);
  double clust_min = clustering[0];
  double clust_max = clustering[1];
  double clust_mean = clustering[2];
  double clust_stdev = clustering[3];
  double clust_log_kurtosis = clustering[4];
  double clust_log_abs_skew = clustering[5];
  double clust_skew_positive = clustering[6];
  double clust_const = clustering[7];
  double clust_time = GetTime(clust_start);

  // Degree metrics
  struct timeval degree_start;
  gettimeofday(&degree_start, 0);
  std::vector<double> degree_metrics;
  GetDegreeData(&degree_metrics);
  double deg_min = degree_metrics[0];
  double deg_max = degree_metrics[1];
  double deg_mean = degree_metrics[2];
  double deg_stdev = degree_metrics[3];
  double deg_log_kurtosis = degree_metrics[4];
  double deg_log_abs_skew = degree_metrics[5];
  double deg_skew_positive = degree_metrics[6];
  double deg_const = degree_metrics[7];
  double degree_time = GetTime(degree_start);
  
  // Compute percentage of edges that have positive weight
  struct timeval ppos_start;
  gettimeofday(&ppos_start, 0);
  double percent_pos = GetPercentPos();
  double ppos_time = GetTime(ppos_start);
  
  // Metrics for weights of graph
  struct timeval weight_start;
  gettimeofday(&weight_start, 0);
  std::vector<double> weights;
  GetWeightData(&weights);
  double weight_min = weights[0];
  double weight_max = weights[1];
  double weight_mean = weights[2];
  double weight_stdev = weights[3];
  double weight_log_kurtosis = weights[4];
  double weight_log_abs_skew = weights[5];
  double weight_skew_positive = weights[6];
  double weight_const = weights[7];
  double weight_time = GetTime(weight_start);
  
  // First two eigenvalues of graph laplacian
  struct timeval ev_start;
  gettimeofday(&ev_start, 0);
  std::pair<double,double> lap_evs = GetLaplacianTopEVs();
  double norm_ev1 = lap_evs.first / avg_degree;
  double norm_ev2 = lap_evs.second / avg_degree;
  double ev_ratio = lap_evs.first / lap_evs.second;
  double log_norm_ev1 = log(std::min(norm_ev1, 1e10));
  double log_norm_ev2 = log(std::min(norm_ev2, 1e10));
  double log_ev_ratio = log(ev_ratio);
  double ev_time  = GetTime(ev_start);

  // Approximation of chromatic number of the graph
  struct timeval chromatic_start;
  gettimeofday(&chromatic_start, 0);
  double chromatic = GetChromaticNumber();
  double chromatic_time = GetTime(chromatic_start);
  
  // Approximation of the diameter, radius, and eccentricity stats of graph
  struct timeval disconnected_start;
  gettimeofday(&disconnected_start, 0);
  double disconnected = Disconnected();
  double disconnected_time = GetTime(disconnected_start);
  
  struct timeval assortativity_start;
  gettimeofday(&assortativity_start, 0);
  double assortativity = DegreeAssortativity();
  double assortativity_time = GetTime(assortativity_start);
  
  struct timeval avg_neighbor_deg_start;
  gettimeofday(&avg_neighbor_deg_start, 0);
  std::vector<double> avg_neighbor_deg;
  AverageNeighborDegree(&avg_neighbor_deg);
  double avg_neighbor_deg_min = avg_neighbor_deg[0];
  double avg_neighbor_deg_max = avg_neighbor_deg[1];
  double avg_neighbor_deg_mean = avg_neighbor_deg[2];
  double avg_neighbor_deg_stdev = avg_neighbor_deg[3];
  double avg_neighbor_deg_log_kurtosis = avg_neighbor_deg[4];
  double avg_neighbor_deg_log_abs_skew = avg_neighbor_deg[5];
  double avg_neighbor_deg_skew_positive = avg_neighbor_deg[6];
  double avg_neighbor_deg_const = avg_neighbor_deg[7];
  double avg_neighbor_deg_time = GetTime(avg_neighbor_deg_start);
  
  struct timeval avg_deg_conn_start;
  gettimeofday(&avg_deg_conn_start, 0);
  std::vector<double> avg_deg_conn;
  AverageDegreeConnectivity(&avg_deg_conn);
  double avg_deg_conn_min = avg_deg_conn[0];
  double avg_deg_conn_max = avg_deg_conn[1];
  double avg_deg_conn_mean = avg_deg_conn[2];
  double avg_deg_conn_stdev = avg_deg_conn[3];
  double avg_deg_conn_log_kurtosis = avg_deg_conn[4];
  double avg_deg_conn_log_abs_skew = avg_deg_conn[5];
  double avg_deg_conn_skew_positive = avg_deg_conn[6];
  double avg_deg_conn_const = avg_deg_conn[7];
  double avg_deg_conn_time = GetTime(avg_deg_conn_start);
  
  struct timeval core_start;
  gettimeofday(&core_start, 0);
  std::vector<double> cores_decomposition;
  CoresDecomposition(&cores_decomposition);
  double core_min = cores_decomposition[0];
  double core_max = cores_decomposition[1];
  double core_mean = cores_decomposition[2];
  double core_stdev = cores_decomposition[3];
  double core_log_kurtosis = cores_decomposition[4];
  double core_log_abs_skew = cores_decomposition[5];
  double core_skew_positive = cores_decomposition[6];
  double core_const = cores_decomposition[7];
  double core_time = GetTime(core_start);
  
  struct timeval mis_start;
  gettimeofday(&mis_start, 0);
  double mis = MaximalIndependentSet();
  double mis_time = GetTime(mis_start);

  if (metrics) {
    *metrics = {clust_min, clust_max, clust_mean, clust_stdev,
                clust_const, deg_min, deg_max, deg_mean, deg_stdev,
                deg_const, percent_pos, weight_min, weight_max,
                weight_mean, weight_stdev, weight_const, chromatic,
                disconnected, assortativity, avg_neighbor_deg_min,
                avg_neighbor_deg_max, avg_neighbor_deg_mean,
                avg_neighbor_deg_stdev, avg_neighbor_deg_const,
                avg_deg_conn_min, avg_deg_conn_max, avg_deg_conn_mean,
                avg_deg_conn_stdev, avg_deg_conn_const, core_min,
                core_max, core_mean, core_stdev, core_const,
                mis, clust_log_kurtosis, deg_log_kurtosis,
                weight_log_kurtosis, avg_neighbor_deg_log_kurtosis,
                avg_deg_conn_log_kurtosis, core_log_kurtosis,
                clust_log_abs_skew, deg_log_abs_skew, weight_log_abs_skew,
                avg_neighbor_deg_log_abs_skew, avg_deg_conn_log_abs_skew,
                core_log_abs_skew, clust_skew_positive, deg_skew_positive,
                weight_skew_positive, avg_neighbor_deg_skew_positive,
                avg_deg_conn_skew_positive, core_skew_positive, log_n,
                log_m, log_norm_ev1, log_norm_ev2, log_ev_ratio};
  }
  if (runtimes) {
    *runtimes = {clust_time, degree_time, ppos_time, weight_time, ev_time,
                 chromatic_time, disconnected_time, assortativity_time,
                 avg_neighbor_deg_time, avg_deg_conn_time, core_time,
                 mis_time};
  }
}

void GraphMetrics::AllMetricNames(std::vector<std::string>* names) {
  if (!names) {
    return;
  }
  *names = {"clust_min", "clust_max", "clust_mean", "clust_stdev",
            "clust_const", "deg_min", "deg_max", "deg_mean", "deg_stdev",
            "deg_const", "percent_pos", "weight_min", "weight_max",
            "weight_mean", "weight_stdev", "weight_const", "chromatic",
            "disconnected", "assortativity", "avg_neighbor_deg_min",
            "avg_neighbor_deg_max", "avg_neighbor_deg_mean",
            "avg_neighbor_deg_stdev", "avg_neighbor_deg_const",
            "avg_deg_conn_min", "avg_deg_conn_max", "avg_deg_conn_mean",
            "avg_deg_conn_stdev", "avg_deg_conn_const", "core_min",
            "core_max", "core_mean", "core_stdev", "core_const",
            "mis", "clust_log_kurtosis", "deg_log_kurtosis",
            "weight_log_kurtosis", "avg_neighbor_deg_log_kurtosis",
            "avg_deg_conn_log_kurtosis", "core_log_kurtosis",
            "clust_log_abs_skew", "deg_log_abs_skew", "weight_log_abs_skew",
            "avg_neighbor_deg_log_abs_skew", "avg_deg_conn_log_abs_skew",
            "core_log_abs_skew", "clust_skew_positive", "deg_skew_positive",
            "weight_skew_positive", "avg_neighbor_deg_skew_positive",
            "avg_deg_conn_skew_positive", "core_skew_positive", "log_n",
            "log_m", "log_norm_ev1", "log_norm_ev2", "log_ev_ratio"};
}

void GraphMetrics::AllRuntimeTypes(std::vector<std::string>* names){
  if (!names) {
    return;
  }
  *names = {"clust_time", "degree_time", "ppos_time", "weight_time", "ev_time",
            "chromatic_time", "disconnected_time", "assortativity_time",
            "avg_neighbor_deg_time", "avg_deg_conn_time", "core_time",
            "mis_time"};
}

void GraphMetrics::GetClusteringData(std::vector<double>* output) {
  // If the graph is complete, then each node clustering coefficient 1
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    std::vector<double> coefs(2, 1.0);  // All clustering coefficient 1
    GetSummary(coefs, output);
    return;
  }

  // Shuffle the nodes and select the numbers we'll look at
  std::vector<int> nodes;  // Shuffled nodes
  for (int i=0; i < n; ++i) {
    nodes.push_back(i);
  }
  std::random_shuffle(nodes.begin(), nodes.end());
  int num_try = std::min<int>(3 * ((int)log(n) + 1), n);
  
  // Compute the clustering coefficient for our selected nodes
  std::vector<double> coefs;  // Custering coefficient of nodes we selected
  for (int idx=0; idx < num_try; ++idx) {
    int i = nodes[idx];  // Node we're computing for
    int num_neighbor = 0;

    // First pass --- build an array of all neighbors
    std::vector<bool> neighbors(n, false);
    for (auto iter = mi_.get_edges_begin(i); iter != mi_.get_edges_end(i);
	 ++iter) {
      neighbors[iter->first] = true;
      ++num_neighbor;
    }

    if (num_neighbor < 2) {
      coefs.push_back(0.0);
    } else {
      // Second pass --- determine triangles (i, j, k) with j < k
      int num_triangle = 0;
      for (auto iter = mi_.get_edges_begin(i); iter != mi_.get_edges_end(i);
           ++iter) {
        for (auto iter2 = mi_.get_edges_begin(iter->first);
             iter2 != mi_.get_edges_end(iter->first); ++iter2) {
          if (iter2->first > iter->first && neighbors[iter2->first]) {
            ++num_triangle;
          }
        }
      }
      coefs.push_back(2.0 * num_triangle / num_neighbor / (num_neighbor - 1));
    }
  }

  // Return the statistics for the selected nodes' clustering coefficients
  GetSummary(coefs, output);
}

void GraphMetrics::GetDegreeData(std::vector<double>* output) {
  int n = mi_.get_size();
  double n_dbl = n - 1.0;  // Normalize by the maximum possible degree of n-1
  std::vector<double> norm_degrees;
  for (int i=0; i < n; ++i) {
    norm_degrees.push_back(mi_.get_vertex_degree(i) / n_dbl);
  }
  GetSummary(norm_degrees, output);
}

double GraphMetrics::GetPercentPos() {
  // Returns the percentage of weights in the graph that are positive.
  int pos_count = 0;
  for (auto iter  = mi_.get_all_edges_begin(); 
            iter != mi_.get_all_edges_end(); iter++) {
    if (iter->second > 0.0) pos_count++;
  }
  return ((double)pos_count)/mi_.get_edge_count();
}


void GraphMetrics::GetWeightData(std::vector<double>* output) {
  std::vector<double> weights;
  double max_abs = 0.0;  // Largest absolute value of a weight
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    double w = iter->second;
    max_abs = std::max<double>(max_abs, fabs(w));
    weights.push_back(w);
  }
  for (auto iter=weights.begin(); iter != weights.end(); ++iter) {
    *iter /= max_abs;
  }
  GetSummary(weights, output);
}


void GraphMetrics::GetSummary(const std::vector<double>& data,
                              std::vector<double>* output) {
  if (!output) {
    return;
  }
  output->clear();

  // Given a vector of data, return a vector of length 4 containing, 
  // in order, [min{data}, max{data}, mean{data}, std{data}]
  double min_x =  std::numeric_limits<double>::max();
  double max_x = -std::numeric_limits<double>::max();
  double sum_x  = 0.0;
  // Pass 1: Compute min, max, and mean of values
  for (auto it = data.begin(); it != data.end(); it++) {
    double d = *it;
    if (d < min_x) min_x = d;
    if (d > max_x) max_x = d;
    sum_x += d;
  }

  // If the data are constant, then exit here, reporting 0 for stdev, skewness,
  // and excess kurtosis. Though skewness and excess kurtosis would actually be
  // infinite in this case, the constant indicator will help us control for the
  // constant case.
  if (min_x == max_x) {
    output->push_back(min_x);  // min
    output->push_back(min_x);  // max
    output->push_back(min_x);  // mean
    output->push_back(0.0);  // stdev
    output->push_back(0.0);  // log_kurtosis
    output->push_back(0.0);  // log_abs_skew
    output->push_back(1.0);  // skew_positive
    output->push_back(1.0);  // constant
    return;
  }

  // Compute second, third and fourth moments of values. While we could have done
  // this in a single iteration above, we'll use a second iteration here for
  // numerical stability.
  double mean = sum_x / data.size();
  double sum_dev2 = 0.0;  // \sum_i (x_i-\bar x)^2
  double sum_dev3 = 0.0;  // \sum_i (x_i-\bar x)^3
  double sum_dev4 = 0.0;  // \sum_i (x_i-\bar x)^4
  for (auto it = data.begin(); it != data.end(); ++it) {
    double dev = *it - mean;
    double dev2 = dev * dev;
    sum_dev2 += dev2;
    sum_dev3 += dev2 * dev;
    sum_dev4 += dev2 * dev2;
  }

  // Compute standard deviation, skewness, and excess kurtosis (plus 3)
  double var         = sum_dev2 / data.size();
  double stdev       = sqrt(var);
  double skewness    = sum_dev3 / data.size() / (stdev*var);
  double ex_kurtosis = sum_dev4 / data.size() / (var*var);

  output->push_back(min_x);
  output->push_back(max_x);
  output->push_back(mean);
  output->push_back(stdev);
  output->push_back(log(ex_kurtosis + 1.0));
  output->push_back(log(fabs(skewness) + 1.0));
  output->push_back(skewness >= 0 ? 1 : 0);
  output->push_back(0.0);  // Are the data constant?
}

double GraphMetrics::Normalize(std::vector<double>* x) {
  // Normalize a vector and return its (original) 2-norm
  int n = x->size();
  double norm = 0.0;
  for (int i=0; i < n; ++i) {
    norm += (*x)[i] * (*x)[i];
  }
  norm = sqrt(norm);
  for (int i=0; i < n; ++i) {
    (*x)[i] /= norm;
  }
  return norm;
}

double GraphMetrics::GetEigenpair(const std::vector<double>& D,
                                  std::vector<double>* x,
                                  std::vector<double>* orthog,
                                  int maxIter, double relDiffLim) {
  // Perform the iterative process described in the comments of
  // GraphMetrics::GetLaplacianTopEVs to compute either the dominant
  // eigenvalue (when orthog is NULL) or the second eigenvalue (when orthog is
  // a pointer to the normalized dominant eigenvector).

  int n = x->size();
  std::vector<double> xp(n);
  double eigenvalue = 0.0;
  double last_eigenvalue = 0.0;

  // Normalize the initial vector
  Normalize(x);

  for (int iter = 0; iter < maxIter; iter++) {
    // Calculate next x' by iterating with graph laplacian and re-orthogonalizing
    xp = D;
    for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
         ++iter) {
      xp[iter->first.first] -= iter->second * (*x)[iter->first.second];
      xp[iter->first.second] -= iter->second * (*x)[iter->first.first];
    }
    if (orthog) {
      double dotProd = 0.0;
      for (int i=0; i < n; ++i) {
        dotProd += xp[i] * (*orthog)[i];
      }
      for (int i=0; i < n; ++i) {
        xp[i] -= dotProd * (*orthog)[i];
      }
    }

    // Calculate norm of x' (this is the eigenvalue) and copy over to x
    last_eigenvalue = eigenvalue;
    eigenvalue = Normalize(&xp);
    *x = xp;

    // Determine if termination criterion reached
    if (iter > 0 &&
        fabs(last_eigenvalue-eigenvalue)/fabs(last_eigenvalue) <= relDiffLim) {
      break;
    }
  }
  return eigenvalue;
}

std::pair<double,double> GraphMetrics::GetLaplacianTopEVs() {
  // The weighted Laplacian matrix of a graph is the matrix L = D - W
  // where W is the weighted adjacency matrix (i.e. W_i,j = weight on i->j) 
  // and D is the (diagonal) degree matrix, D_i,i = sum_j W_i,j
  //
  // Power method of finding dominant eigenvalue:
  // To find the dominant eigenvalue of the matrix A, we initialize a random
  // vector x. We then calculate a vector x' = Ax. We set x = x'/||x'|| and
  // repeat the matrix-vector operation. As we continue to do this we
  // find that ||x'|| converges to the dominant eigenvalue.
  //
  // To find the second eigenvalue, we'll again start with a random vector and
  // calculate x' = Ax. We then orthogonalize with x' = x' - <x', x1>x1, where
  // x1 is the normalized dominant eigenvalue. Finally we set x = x'/||x'|| and
  // repeat the process. As we continue to do this, ||x'|| converges to the
  // second eigenvalue.
  //
  // We apply the power method to the weighted Laplacian, i.e.
  // Lx = (D - W)x = Dx - Wx
  // We'll compute D once, but store it as a vector (i.e. we'll never form
  // L explicitly).
  int n = mi_.get_size();

  // Precalculate D
  std::vector<double> D(n, 0.0);
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    D[iter->first.first] += iter->second;
    D[iter->first.second] += iter->second;
  }

  // Setup memory for x vectors and initalize to random
  std::vector<double> x1(n);  // First eigenvector
  std::vector<double> x2(n);  // Second eigenvector
  for (int i = 0; i < n; i++) {
    x1[i] = Random::RandDouble();
    x2[i] = Random::RandDouble();
  }

  // Compute first two eigenvalues
  double e1 = GetEigenpair(D, &x1, NULL, 10, 1e-6);
  double e2 = GetEigenpair(D, &x2, &x1, 10, 1e-6);
  return std::pair<double,double>(e1, e2);
}


double GraphMetrics::GetChromaticNumber() {
  // Use a greedy heurisic (Welsh-Powell) to color the graph
  // http://en.wikipedia.org/wiki/Greedy_coloring
  // This was helpful for implementing
  // http://mrsleblancsmath.pbworks.com/w/file/fetch/46119304/vertex%20coloring%20algorithm.pdf

  int n = mi_.get_size();

  // On a complete graph the chromatic number is just the number of nodes
  if (mi_.get_edge_count() == n*(n-1)/2) {
    return 1.0;  // Return value is normalized by number of nodes
  }

  // Descending order by degree
  std::vector<std::pair<int,int>> degs;
  for (int i = 0; i < n; i++) {
    degs.push_back(std::pair<int,int>(-mi_.get_vertex_degree(i),i));
  }
  std::sort(degs.begin(), degs.end());

  std::vector<bool> colored(n, false);  // Has a node been colored?
  int num_rem = n;  // How many nodes still need to be colored?
  int chromatic_num = 0;  // Number of colors needed in total
  while (num_rem > 0) {
    ++chromatic_num;

    // Is each node adjacent to a node colored this iteration?
    std::vector<bool> adjacent(n, false);

    // For every node (in order of degree), try to color
    for (int ind=0; ind < n; ++ind) {
      int i = degs[ind].second;
      if (colored[i] || adjacent[i]) {
        continue;  // Already colored or adjacent to a node with this iter's col.
      }

      // Color this node and mark all its neighbors as being adjacent
      colored[i] = true;
      --num_rem;
      for (auto iter=mi_.get_edges_begin(i); iter != mi_.get_edges_end(i);
           ++iter) {
        adjacent[iter->first] = true;
      }
    }
  }
  
  return ((double)chromatic_num) / n;  // Return normalized chromatic num
}

bool GraphMetrics::Disconnected() {
  // Determine if the graph is disconnected by performing a DFS from node 0.
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    // Complete graphs are not disconnected
    return false;
  }

  std::vector<bool> added(n, false);  // Have we reached this node in DFS?
  std::vector<int> toVisit;
  toVisit.push_back(0);  // Start at node 0
  while (toVisit.size() > 0) {
    int node = toVisit.back();
    toVisit.pop_back();
    added[node] = true;
    for (auto iter = mi_.get_edges_begin(node); iter != mi_.get_edges_end(node);
         ++iter) {
      if (!added[iter->first]) {
        toVisit.push_back(iter->first);
      }
    }
  }

  // If we've visited all the nodes, then we're connected. Otherwise, we're
  // disconnected.
  for (int i=0; i < n; ++i) {
    if (!added[i]) {
      return true;
    }
  }
  return false;
}

double GraphMetrics::DegreeAssortativity() {
  // For each edge (i, j) and associated node degrees d_i and d_j, append
  // d_i and d_j to a vector x in that order and append d_j and d_i to a vector
  // y in that order. Return the Pearson's correlation coefficient between the
  // vectors. For numerical stability, we'll store the vectors and do a second
  // pass to compute the correlation. Note that assortativity is typically stated
  // as being computed with the degree minus one, but this yields an identical
  // correlation.
  
  // In a complete graph, just return 0 (degrees are constant on all edges, and
  // we return 0 in this case below).
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    return 0.0;
  }

  // Grab degrees of each node
  std::vector<int> degrees;
  for (int i=0; i < mi_.get_size(); ++i) {
    degrees.push_back(mi_.get_vertex_degree(i));
  }

  // Pass 1: Build x and y, retaining the sum
  std::vector<int> x;
  std::vector<int> y;
  int sum = 0;  // Sum of elements in x (and by symmetry y)
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    x.push_back(degrees[i]);
    x.push_back(degrees[j]);
    y.push_back(degrees[j]);
    y.push_back(degrees[i]);
    sum += degrees[i] + degrees[j];
  }
  double mean = ((double)sum) / x.size();

  // Pass 2: Compute correlation coefficient
  double sum_dxdx = 0.0;  // Sum of (x-mean)^2; equivalent to sum of (y-mean)^2
  double sum_dxdy = 0.0;  // Sum of (x-mean)(y-mean)
  for (int i=0; i < x.size(); ++i) {
    double dx = x[i] - mean;
    double dy = y[i] - mean;
    sum_dxdx += dx * dx;
    sum_dxdy += dx * dy;
  }

  // Return the correlation coefficient (0 if degrees are constant)
  if (sum_dxdx > 0.0) {
    return sum_dxdy / sum_dxdx;
  } else {
    return 0.0;
  }
}

void GraphMetrics::AverageNeighborDegree(std::vector<double>* output) {
  // Compute the average neighbor degree of each node (0 if no neighbors;
  // normalized by the maximum possible neighbor count) and call GetSummary on
  // that.

  // In complete graphs, constant with normalized value 1
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    std::vector<double> norm_degs(2, 1.0);
    GetSummary(norm_degs, output);
    return;
  }
  
  // Get the normalized degree of each node
  double n_dbl = n - 1.0;  // Normalize by the maximum possible degree of n-1
  std::vector<double> norm_degrees;
  for (int i=0; i < n; ++i) {
    norm_degrees.push_back(mi_.get_vertex_degree(i) / n_dbl);
  }

  // For each node with one or more neighbors, compute their average degree
  std::vector<double> neighbor_degree_sum(n, 0.0);  // Sum of degs for each node
  std::vector<int> neighbor_count(n, 0);  // Number of neighbors
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    int i = iter->first.first;
    int j = iter->first.second;
    neighbor_degree_sum[i] += norm_degrees[j];
    neighbor_degree_sum[j] += norm_degrees[i];
    ++(neighbor_count[i]);
    ++(neighbor_count[j]);
  }

  // Compute average normalized degree of neighbors (ignore nodes with no
  // neighbors)
  std::vector<double> avg_neighbor_degree;
  for (int i=0; i < n; ++i) {
    if (neighbor_count[i] > 0) {
      avg_neighbor_degree.push_back(neighbor_degree_sum[i] / neighbor_count[i]);
    }
  }
  
  // Return summary statistics for the distribution (constant 0 if no edges)
  if (avg_neighbor_degree.size() == 0) {
    avg_neighbor_degree.push_back(0.0);
  }
  GetSummary(avg_neighbor_degree, output);
}

void GraphMetrics::AverageDegreeConnectivity(std::vector<double>* output) {
  // Compute the average neighbor degree of nodes with each degree count
  // (normalized by the maximum possible neighbor count) and call GetSummary on
  // that.

  // In complete graphs, constant with normalized value 1
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    std::vector<double> norm_degs(2, 1.0);
    GetSummary(norm_degs, output);
    return;
  }
  
  // Get the degree of each node
  double n_dbl = n - 1.0;  // Normalize by the maximum possible degree of n-1
  std::vector<int> degrees;
  for (int i=0; i < n; ++i) {
    degrees.push_back(mi_.get_vertex_degree(i));
  }

  // For each node with one or more neighbors, compute their average degree
  std::vector<double> neighbor_degree_sum(n, 0.0);  // Sum of degs for each deg
  std::vector<int> neighbor_count(n, 0);  // Number of neighbors
  for (auto iter = mi_.get_all_edges_begin(); iter != mi_.get_all_edges_end();
       ++iter) {
    int deg_i = degrees[iter->first.first];
    int deg_j = degrees[iter->first.second];
    neighbor_degree_sum[deg_i] += deg_j / n_dbl;
    neighbor_degree_sum[deg_j] += deg_i / n_dbl;
    ++(neighbor_count[deg_i]);
    ++(neighbor_count[deg_j]);
  }

  // Compute average normalized degree of neighbors (ignore degrees that don't
  // show up in our graph)
  std::vector<double> avg_neighbor_degree;
  for (int deg=1; deg < n; ++deg) {
    if (neighbor_count[deg] > 0) {
      avg_neighbor_degree.push_back(neighbor_degree_sum[deg] /
                                    neighbor_count[deg]);
    }
  }
  
  // Return summary statistics for the distribution (constant 0 if no edges)
  if (avg_neighbor_degree.size() == 0) {
    avg_neighbor_degree.push_back(0.0);
  }
  GetSummary(avg_neighbor_degree, output);
}

void GraphMetrics::CoresDecomposition(std::vector<double>* output) {
  // Compute the core number for each node in the graph, following
  // http://arxiv.org/pdf/cs/0310049v1.pdf. Normalize by the maximum core number,
  // n-1. Return GetSummary called on this vector of numbers.

  // In complete graphs, constant with normalized value 1
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    std::vector<double> norm_degs(2, 1.0);
    GetSummary(norm_degs, output);
    return;
  }

  // Get the degree of each node, as well as a maximum degree
  std::vector<int> degrees;
  int max_degree = 0;
  for (int i=0; i < n; ++i) {
    degrees.push_back(mi_.get_vertex_degree(i));
    max_degree = std::max<int>(max_degree, degrees[i]);
  }
  
  // Construct the following matrices (lines 13-27 of Alg. 1 in linked paper)
  // bin[0:max_degree] -- bin[i] is location of first element of degree >= i
  // vert[0:n] -- sorted list of vertices by degree
  // pos[0:n] -- pos[i] is position of node i in vert
  std::vector<int> bin(max_degree+1, 0);
  std::vector<int> vert(n);
  std::vector<int> pos(n);
  for (int i=0; i < n; ++i) {
    ++(bin[degrees[i]]);
  }
  int start = 0;
  for (int deg=0; deg <= max_degree; ++deg) {
    int num = bin[deg];
    bin[deg] = start;
    start += num;
  }
  for (int i=0; i < n; ++i) {
    pos[i] = bin[degrees[i]];
    vert[pos[i]] = i;
    ++(bin[degrees[i]]);
  }
  for (int deg=max_degree; deg >= 1; --deg) {
    bin[deg] = bin[deg-1];
  }
  bin[0] = 0;

  // Compute the cores of each node
  for (int i=0; i < n; ++i) {
    int v = vert[i];
    for (auto iter = mi_.get_edges_begin(v); iter != mi_.get_edges_end(v); ++iter) {
      int u = iter->first;  // Neighbor of v
      int du = degrees[u];  // Degree of u
      if (du > degrees[v]) {
        // We will decrease the degree of u by 1 and move it up in vert if needed
        int pu = pos[u];  // Position of u in vert
        int pw = bin[du];  // Position where u will be moved
        int w = vert[pw];  // Vertex to be swapped with u
        if (u != w) {
          // Flip u and w
          pos[u] = pw;
          vert[pw] = u;
          pos[w] = pu;
          vert[pu] = w;
        }
        ++(bin[du]);
        --(degrees[u]);
      }
    }
  }

  // degrees now stores the cores; normalize by n-1 and get statistics
  std::vector<double> norm_degrees;
  for (int i=0; i < n; ++i) {
    norm_degrees.push_back(((double)degrees[i]) / (n-1));
  }
  GetSummary(norm_degrees, output);
}

double GraphMetrics::MaximalIndependentSet() {
  // A maximal independent set is a maximal collection of vertices where no two
  // vertices are adjacent. We normalize our result by n, the number of nodes.

  // In a complete graph normalized output is 1/n
  int n = mi_.get_size();
  if (mi_.get_edge_count() == n*(n-1)/2) {
    return 1.0 / n;
  }

  // Sort nodes by degree (increasing) so we add lowest-degree nodes first
  std::vector<std::pair<int, int> > nodes;
  for (int i=0; i < n; ++i) {
    nodes.push_back(std::pair<int, int>(mi_.get_vertex_degree(i), i));
  }
  std::sort(nodes.begin(), nodes.end());

  std::vector<bool> adj(n, false);  // Adjacent to a node in maximal set?
  int misSize = 0;  // Size of current maximal independent set
  for (int idx=0; idx < n; ++idx) {
    // Loop through nodes in increasing order of degree. If a node is not
    // adjacent to something already in the maximal set, add it to maximal set
    // and label its neighbors as being adjacent
    int i = nodes[idx].second;
    if (!adj[i]) {
      ++misSize;
      for (auto iter = mi_.get_edges_begin(i); iter != mi_.get_edges_end(i);
           ++iter) {
        adj[iter->first] = true;
      }
    }
  }

  return ((double)misSize) / n;
}
