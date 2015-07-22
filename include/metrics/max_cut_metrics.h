#ifndef METRICS_MAX_CUT_METRICS_H_
#define METRICS_MAX_CUT_METRICS_H_

#include <string>
#include <vector>
#include "problem/max_cut_instance.h"

class GraphMetrics {
 public:
  GraphMetrics(const MaxCutInstance& mi);

  void AllMetrics(std::vector<double>* metrics, std::vector<double>* runtimes);
  static void AllMetricNames(std::vector<std::string>* names);
  static void AllRuntimeTypes(std::vector<std::string>* names);

  void GetClusteringData(std::vector<double>* output);
  void GetDegreeData(std::vector<double>* output);
  double GetPercentPos();
  void GetWeightData(std::vector<double>* output);

  std::pair<double,double> GetLaplacianTopEVs();

  double GetChromaticNumber();
  bool Disconnected();

  double DegreeAssortativity();
  void AverageNeighborDegree(std::vector<double>* output);
  void AverageDegreeConnectivity(std::vector<double>* output);
  void CoresDecomposition(std::vector<double>* output);
  double MaximalIndependentSet();

 private:
  // Max-Cut instance we're operating on
  const MaxCutInstance& mi_;

  // Helper functions
  static double GetTime(const struct timeval& start);
  static void GetSummary(const std::vector<double>& data,
                         std::vector<double>* output);
  static double Normalize(std::vector<double>* x);
  double GetEigenpair(const std::vector<double>& D,
                      std::vector<double>* x, std::vector<double>* orthog,
                      int maxIter, double relDiffLim);
};

#endif
