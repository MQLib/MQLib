#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <iostream>
#include <string>

#include "heuristics/heuristic_factory.h"
#include "heuristics/maxcut/hyperheuristic.h"
#include "metrics/max_cut_metrics.h"
#include "problem/max_cut_instance.h"
#include "problem/qubo_instance.h"
#include "util/ezOptionParser.h"

void Usage(ez::ezOptionParser& opt) {
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
};

double GetTime(const struct timeval& start) {
  // Runtime in seconds since the passed start time
  struct timeval end;
  gettimeofday(&end, 0);
  return (end.tv_sec - start.tv_sec) + 0.000001 * (end.tv_usec - start.tv_usec);
}

// The runtime limit according to instance size (may be overwritten on command
// line).
double RuntimeLimit(MaxCutInstance *mi, QUBOInstance *qi) {
  int n = mi ? mi->get_size() : (qi->get_size() + 1);
  double runtime_limit = 0.59 * n;
  if (runtime_limit < 120) {
    runtime_limit = 120;
  } else if (runtime_limit > 1200) {
    runtime_limit = 1200;
  }
  return runtime_limit;
}

int main(int argc, const char* argv[]) {
  ez::ezOptionParser opt;

  opt.overview = "MQLib: Library of Max-Cut and QUBO heuristics";
  opt.syntax = "\n# Run Max-Cut or QUBO heuristic\n./bin/MQlib -h heur_code | -hh -fM maxcut_file [-nv] [-ps] [-q | -r runtime_limit] [-s SEED]\n./bin/MQlib -h heur_code | -hh -fQ qubo_file [-nv] [-ps] [-q | -r runtime_limit] [-s SEED]\n\n# Compute metrics for an input file\n./bin/MQlib -fM maxcut_file [-mh] [-m]\n./bin/MQlib -fQ qubo_file [-mh] [-m]\n\n# List the available heuristics.\n./bin/MQlib -l";
  opt.example = "./bin/MQlib -h BURER2002 -fM bin/sampleMaxCut.txt -r 10\n";

  opt.add("",  // Default
	  0,  // Required?
	  1,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Heuristic code.",  // Help description
	  "-h",  // Flag token
	  "--heuristic"
	  );
  
  opt.add("",  // Default
          0,  // Required?
          0,  // Number of args expected
          0,  // Delimiter if expecting multiple args
          "Use the Max-Cut hyper-heuristic",  // Help description
          "-hh",  // Flag token
          "--hyperheuristic"
          );

  opt.add("",  // Default
	  0,  // Required?
	  1,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Filename for Max-Cut problem instance.",  // Help description
	  "-fM",  // Flag token
	  "--fileMaxCut"
	  );

  opt.add("",  // Default
	  0,  // Required?
	  1,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Filename for QUBO problem instance.",  // Help description
	  "-fQ",  // Flag token
	  "--fileQUBO"
	  );

  // Limit seed to range of unsigned short
  ez::ezOptionValidator* vU2 = new ez::ezOptionValidator("u2");
  opt.add("",  // Default
	  0,  // Required?
	  1,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Random seed (range 0 to 65535).",  // Help description
	  "-s",  // Flag token
	  "--seed",
	  vU2
	  );

  opt.add("",  // Default
	  0,  // Required?
	  0,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "No validation (after the run is complete, don't check the solutions to make sure objectives were reported accurately).",  // Help description
	  "-nv",  // Flag token
	  "--no-validation"
	  );

  opt.add("",  // Default
	  0,  // Required?
	  0,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Output problem instance statistics.",  // Help description
	  "-m",  // Flag token
	  "--metrics"
	  );

  opt.add("",  // Default
          0,  // Required?
          0,  // Number of args expected
          0,  // Delimiter if expecting multiple args
          "Output header for problem instance statistics",  // Help description
          "-mh",  // Flag token
          "--metricsHeader"
          );

  opt.add("",  // Default
          0,   // Required?
          0,   // Number of args expected
          0,   // Delimiter if expecting multiple args
          "Print solution to screen",  // Help description
          "-ps",  // Flag token
          "--printSolution"
          );

  opt.add("",  // Default
	  0,  // Required?
	  0,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Fast run (10x faster than full run)",  // Help description
	  "-q",  // Flag token
	  "--quick"
	  );

  opt.add("",  // Default
          0,  // Required?
          0,  // Number of args expected
          0,  // Delimiter if expecting multiple args
          "List all heuristic codes and brief descriptions",  // Help description
          "-l",  // Flag token
          "--listHeuristics"
          );

  // Runtime limit must be a double
  ez::ezOptionValidator* vD = new ez::ezOptionValidator("d");
  opt.add("",  // Default
	  0,  // Required?
	  1,  // Number of args expected
	  0,  // Delimiter if expecting multiple args
	  "Runtime limit (seconds), or iteration count for baseline",  // Help description
	  "-r",
	  "--runtime",
	  vD
	  );

  opt.parse(argc, argv);

  std::vector<std::string> badOptions;
  if (!opt.gotRequired(badOptions)) {
    for (int i=0; i < badOptions.size(); ++i) {
      std::cout << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
    }
    Usage(opt);
    return 1;
  }

  if (!opt.gotExpected(badOptions)) {
    for (int i=0; i < badOptions.size(); ++i) {
      std::cout << "ERROR: Got unexpected number of arguments for option " <<
	badOptions[i] << ".\n\n";
    }
    Usage(opt);
    return 1;
  }

  std::vector<std::string> badArgs;
  if (!opt.gotValid(badOptions, badArgs)) {
    for (int i=0; i < badOptions.size(); ++i) {
      std::cerr << "ERROR: Got invalid argument \"" << badArgs[i] <<
	"\" for option " << badOptions[i] << ".\n\n";
    }
    Usage(opt);
    return 1;
  }

  // Check if any of the options for a heuristic run are set
  bool heurSet = opt.isSet("-h") || opt.isSet("-hh") || opt.isSet("-nv") ||
    opt.isSet("-ps") || opt.isSet("-q") || opt.isSet("-r") || opt.isSet("-s");
  bool metricSet = opt.isSet("-m") || opt.isSet("-mh");
  bool listSet = opt.isSet("-l");
  int numSet = ((int)heurSet) + ((int)metricSet) + ((int)listSet);
  if (numSet != 1) {
    std::cout << "ERROR: Invalid usage." << std::endl;
    Usage(opt);
    return 1;
  }

  /*** Handle listSet case ***/
  if (listSet) {
    HeuristicFactory factory;
    factory.PrintHeuristicCodes();
    return 0;
  }

  // -q and -r should not be used together
  if (opt.isSet("-q") && opt.isSet("-r")) {
    std::cout << "ERROR: -q and -r should not be used together" << std::endl;
    Usage(opt);
    return 1;
  }

  // Exactly one of -h and -hh is required when running a heuristic
  if (heurSet && !opt.isSet("-h") && !opt.isSet("-hh")) {
    std::cout << "ERROR: You must provide a heuristic code with -h or -hh " <<
      "when running a heuristic." << std::endl;
    Usage(opt);
    return 1;
  }
  if (opt.isSet("-h") && opt.isSet("-hh")) {
    std::cout << "ERROR: -h and -hh cannot both be used" << std::endl;
    Usage(opt);
    return 1;
  }

  // Exactly one of -fM and -fQ required unless we are running with -mh alone
  // in which case you can run with neither set.
  if (!opt.isSet("-fM") && !opt.isSet("-fQ") &&
      (!metricSet || opt.isSet("-m"))) {
    std::cout << "ERROR: Either -fM or -fQ must be used" << std::endl;
    Usage(opt);
    return 1;
  }
  if (opt.isSet("-fM") && opt.isSet("-fQ")) {
    std::cout << "ERROR: -fM and -fQ cannot both be provided" << std::endl;
    Usage(opt);
    return 1;
  }

  // Build the problem instance
  std::string filename;
  MaxCutInstance* mi = 0;
  QUBOInstance* qi = 0;
  if (opt.isSet("-fM")) {
    opt.get("-fM")->getString(filename);
    mi = new MaxCutInstance(filename);
  } else if (opt.isSet("-fQ")) {
    opt.get("-fQ")->getString(filename);
    qi = new QUBOInstance(filename);
  }

  if (heurSet) {
    /************* Handle heuristicSet case ****************/
    // Set the seed if it's provided; otherwise set to current time
    int seed = time(0);
    if (opt.isSet("-s")) {
      opt.get("-s")->getInt(seed);
    }
    srand(seed);
    
    // Compute the runtime limit
    double runtime_limit = RuntimeLimit(mi, qi);
    if (opt.isSet("-q")) {
      runtime_limit *= 0.1;
    } else if (opt.isSet("-r")) {
      opt.get("-r")->getDouble(runtime_limit);
      if (runtime_limit < 0.0) {
        std::cout << "Illegal runtime limit: " << runtime_limit << " seconds" <<
          std::endl;
        return 1;
      }
    }
    
    // Run the heuristic
    bool validation = !opt.isSet("-nv");
    MaxCutHeuristic *mh = NULL;
    QUBOHeuristic *qh = NULL;
    Heuristic* heuristic = NULL;
    std::string heuristic_code;
    if (opt.isSet("-h")) {
      // Run a specified heuristic
      opt.get("-h")->getString(heuristic_code);
      HeuristicFactory factory;
      if (factory.ValidMaxCutHeuristicCode(heuristic_code)) {
        if (!mi) {
          mi = new MaxCutInstance(*qi);
        }
        mh = factory.RunMaxCutHeuristic(heuristic_code, *mi, runtime_limit,
                                        validation, NULL);
        heuristic = mh;
      } else if (factory.ValidQUBOHeuristicCode(heuristic_code)) {
        if (!qi) {
          qi = new QUBOInstance(*mi);
        }
        qh = factory.RunQUBOHeuristic(heuristic_code, *qi, runtime_limit,
                                      validation, NULL);
        heuristic = qh;
      } else {
        std::cout << "Illegal heuristic code " << heuristic_code << std::endl;
        std::cout << std::endl;
        factory.PrintHeuristicCodes();
        return 1;
      }
    } else if (opt.isSet("-hh")) {
      // Run the Max-Cut hyperheuristic
      if (!mi) {
        mi = new MaxCutInstance(*qi);
      }
      std::string selected;
      mh = new MaxCutHyperheuristic(*mi, runtime_limit, validation, NULL, seed,
                                    &selected);
      heuristic = mh;
      heuristic_code = "HH_" + selected;
    }


    // Freeze true run time before validation
    const double final_runtime = heuristic->Runtime();
    
    // Validate solutions and output data if valid
    if (validation && !heuristic->IsHistoryValid()) {
      std::cout << "Error: heuristic history was invalid!" << std::endl;
    } else {
      // Common output whenever we're running heuristics
      std::cout << runtime_limit << "," << heuristic_code << ",\"" << filename <<
        "\"," << std::setprecision(15) << heuristic->get_best() << "," <<
        final_runtime << "," << heuristic->History() << std::endl;
    }
    
    // Print out the final solution if requested
    if (opt.isSet("-ps")) {
      std::cout << std::endl << "Solution:" << std::endl;
      if (mh && opt.isSet("-fM")) {
        const MaxCutSimpleSolution& sol = mh->get_best_solution();
        sol.PrintSolution();
      } else if (mh && opt.isSet("-fQ")) {
        QUBOSimpleSolution sol(mh->get_best_solution(), *qi, NULL);
        sol.PrintSolution();
      } else if (qh && opt.isSet("-fM")) {
        MaxCutSimpleSolution sol(qh->get_best_solution(), *mi, NULL);
        sol.PrintSolution();
      } else if (qh && opt.isSet("-fQ")) {
        const QUBOSimpleSolution& sol = qh->get_best_solution();
        sol.PrintSolution();
      }
    }

    if (mh) {
      delete mh;
      mh = NULL;
    }
    if (qh) {
      delete qh;
      qh = NULL;
    }
  } else {
    /************* Handle metricSet case ****************/
    if (opt.isSet("-mh")) {
      // Column order is [metrics], [runtimes]
      std::vector<std::string> metric_names;
      GraphMetrics::AllMetricNames(&metric_names);
      for (int i=0; i < metric_names.size(); ++i) {
        if (i != 0) {
          std::cout << ",";
        }
        std::cout << metric_names[i];
      }
      std::vector<std::string> runtime_names;
      GraphMetrics::AllRuntimeTypes(&runtime_names);
      for (int i=0; i < runtime_names.size(); ++i) {
        std::cout << "," << runtime_names[i];
      }
      std::cout << std::endl;
    }
    if (opt.isSet("-m")) {
      // Output metrics associated with this problem instance
      // Convert QUBO instance to Max-Cut instance if needed
      if (!mi) {
        mi = new MaxCutInstance(*qi);
      }

      // Column order is [metrics], [runtimes]
      std::vector<double> metrics;
      std::vector<double> runtimes;
      GraphMetrics gm(*mi);
      gm.AllMetrics(&metrics, &runtimes);
      for (int i=0; i < metrics.size(); ++i) {
        if (i != 0) {
          std::cout << ",";
        }
        std::cout << metrics[i];
      }
      for (int i=0; i < runtimes.size(); ++i) {
        std::cout << "," << runtimes[i];
      }
      std::cout << std::endl;
    }
  }

  // Clean up memory allocated for instances
  if (mi) {
    delete mi;
    mi = 0;
  }
  if (qi) {
    delete qi;
    qi = 0;
  }

  return 0;
}
