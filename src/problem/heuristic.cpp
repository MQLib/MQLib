#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include "problem/heuristic.h"

Heuristic::Heuristic(double runtime_limit, bool validation) :
  validation_(validation),
  best_(0.0),
  runtime_limit_(runtime_limit) {
  gettimeofday(&start_time_, 0);
}

double Heuristic::Runtime() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  double secs = (tv.tv_sec - start_time_.tv_sec) + 0.000001 *
    (tv.tv_usec - start_time_.tv_usec);
  return secs;
}

std::string Heuristic::History() {
  std::stringstream out_str;
  out_str << std::setprecision(15) << "[";
  for (int i = 0; i < past_solution_values_.size(); i++) {
    if (i > 0) out_str << ";";
    out_str << past_solution_values_[i] << ":" << past_solution_times_[i];
  }
  out_str << "]";
  return out_str.str();
}
