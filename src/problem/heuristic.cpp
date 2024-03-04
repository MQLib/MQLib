#include <iomanip>
#include <iostream>
#include "mqlib/problem/heuristic.h"

namespace mqlib {

    Heuristic::Heuristic(double runtime_limit, bool validation) :
            validation_(validation),
            best_(0.0),
            runtime_limit_(runtime_limit) {
        gettimeofday(&start_time_, nullptr);
    }

    double Heuristic::Runtime() const {
        struct timeval tv;
        gettimeofday(&tv, nullptr);
        double secs = static_cast<double>(tv.tv_sec - start_time_.tv_sec) + 0.000001 *
                                                                            static_cast<double>(tv.tv_usec -
                                                                                                start_time_.tv_usec);
        return secs;
    }

    std::string Heuristic::History() {
        std::stringstream out_str;
        out_str << std::setprecision(15) << "[";
        for (uint64_t i = 0; i < past_solution_values_.size(); i++) {
            if (i > 0) out_str << ";";
            out_str << past_solution_values_[i] << ":" << past_solution_times_[i];
        }
        out_str << "]";
        return out_str.str();
    }

    const std::vector<double> &Heuristic::get_past_solution_values() const {
        return past_solution_values_;
    }

    const std::vector<double> &Heuristic::get_past_solution_times() const {
        return past_solution_times_;
    }

}
