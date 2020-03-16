# Using the MQLib

This guide describes how to use the MQLib to run heuristics and compute problem instance metrics. We first describe how to run heuristics and compute problem-instance metrics through the command-line interface, next describe the file format used to specify problem instances, and finally describe how to run a heuristic by linking to the MQLib.

## Running MQLib From the Command Line

After building MQLib (see the [MQLib guide](../README.md) for details on this step), an executable providing a command-line interface to the library is available at `bin/MQLib`. This executable can be used for two main purposes: running a heuristic and outputting metrics for a Max-Cut problem instance.

### Running a heuristic

You can run a heuristic by specifying the heuristic to be run and the instance to be run on (either a Max-Cut instance or a QUBO instance). If a Max-Cut heuristic is specified for a QUBO instance, then the instance will be reduced to a Max-Cut instance, the heuristic will be used, and then the solutions will be converted to solutions for the original QUBO instance. Similarly, if a QUBO heuristic is specified for a Max-Cut instance, the instance will be reduced to a QUBO instance, the heuristic will be run, and solutions will be converted back.

As an example, if you wanted to run Max-Cut heuristic `BURER2002` on the sample QUBO instance provided with the repository, [bin/sampleQUBO.txt](sampleQUBO.txt), for 10 seconds, outputting the final solution, you could run `bin/MQLib -fQ bin/sampleQUBO.txt -h BURER2002 -r 10 -ps` from the main MQLib folder. This creates output of the following format:

```
10,BURER2002,"bin/sampleQUBO.txt",6,10.00002,[0:0.000418;6:0.000797]

Solution:
1 0 1
```

In order, the elements of the first line are:

1. The runtime limit
2. The code of the heuristic run
3. The name of the instance run on
4. The best objective value of any solution found by the heuristic
5. The total time the heuristic ran before terminating
6. The history of all new best solutions returned by the heuristic.

In this case, the heuristic began running after 0.000418 seconds of startup time, found a new best solution with objective value 6 at time 0.000797 seconds, and then continued (unsuccessfully) searching for better solutions for the remaining computational budget.

A number of command-line options are available to customize runs:
* `-h` / `-hh`: Specifies the heuristic to be run. Flag `-h` specifies the name of a heuristic to run. The list of all available heuristics (along with a brief description of each) is available by running `bin/MQLib -l` from the main MQLib folder. Flag `-hh` specifies that the hyper-heuristic should be run. The hyper-heuristic selects the heuristic to be run on the instance based on the instance's properties. The run above could be changed to use the hyperheuristic with `bin/MQLib -fQ bin/sampleQUBO.txt -hh -r 10 -ps`.
* `-fM` / `-fQ`: Specifies the name of a file describing a Max-Cut (QUBO) instance. See later in this README for a description of file formats.
* `-nv`: Turns off the validation check that is performed after the heuristic run is complete. This validation check verifies that each new best solution has an accurate objective value based on its variable values.
* `-ps`: Print the best solution found.
* `-r` / `-q`: If `-r` is specified, then this is the runtime limit, in seconds. If `-r` is omitted, then the runtime limit is set to `0.59*n`, where `n` is the number of nodes in the instance (or the number of QUBO variables, plus one). This runtime limit is then clamped to be no smaller than 120 seconds and no larger than 1200 seconds. If `-q` is specified, then the total runtime is one tenth of this computed runtime limit.
* `-s`: The random number generator seed to be used for the run.

### Compute metrics for a Max-Cut problem instance

You can compute metrics for a Max-Cut problem instance using the `-m` flag. For instance, to compute the metrics associated with the sample Max-Cut instance provided with the repository, [bin/sampleMaxCut.txt](sampleMaxCut.txt), one would run `bin/MQLib -fM bin/sampleMaxCut.txt -m` from the main MQLib folder. The resulting output is the set of all metrics calculated for the problem instance, in csv format. To include a header with the name of each metric, the `-mh` flag can also be provided.

If metrics are instead requested for a QUBO instance (indicated by the `-fQ` flag instead of the `-fM` flag), then the instance is reduced to a Max-Cut instance before the metrics are computed.

## File Format

We use a simple file format, in which any line that begins with a `#` is treated as a comment. The first non-comment line is of the form `n m`, where `n` is the number of nodes in the instance (Max-Cut) or the number of variables in the instance (QUBO). The other number is the number of subsequent "data lines" in the file.

Each subsequent line is a "data line" of the form `a b w`. In a Max-Cut instance, this refers to an edge between nodes `a` and `b` of weight `w` (edges should only be added once, so there should be no corresponding `b a w` line -- this is important, since we do very limited integrity checking on files). In a QUBO instance, this refers to an entry `w = Q_ab` in the input matrix `Q`. Similarly, if there is a line `a b w` is included then the line `b a w` should not be included. In both cases, `a` and `b` are 1-indexed, meaning they vary from 1 to `n`.

An example Max-Cut instance (the one in [bin/sampleMaxCut.txt](sampleMaxCut.txt)) might have three nodes, with nodes A and B being connected by an edge of weight 3 and nodes B and C being connected by an edge of weight -2. This would be encoded as:

```
3 2
1 2 3
2 3 -2
```

An example QUBO instance (the one in [bin/sampleQUBO.txt](sampleQUBO.txt)) might have three variables, with input matrix 

```
[ 5 -6  0
 -6  3 -1
  0 -1  1]
```

This would be encoded as:

```
3 5
1 1 5
2 2 3
3 3 1
1 2 -6
2 3 -1
```

## Linking to MQLib

While building the MQLib, the library `bin/MQLib.a` should be generated. This can be used to link to the MQLib. As the MQLib is released under the MIT license (see the [LICENSE](../LICENSE) file), such linking should not be problematic for most software projects.

### Basic Heuristic Run

To illustrate the basic code needed to run a heuristic and get information about the best solution found, we will implement a program (which we'll assume is stored in a file `linked.cpp`) that runs the Max-Cut heuristic `BURER2002` on the sample Max-Cut problem instance stored in [bin/sampleMaxCut.txt](sampleMaxCut.txt), running for 10 seconds:

```
#include <iostream>
#include "heuristics/maxcut/burer2002.h"

int main(int argc, char** argv) {
  MaxCutInstance mi("bin/sampleMaxCut.txt");
  Burer2002 heur(mi, 10.0, false, NULL);
  const MaxCutSimpleSolution& mcSol = heur.get_best_solution();
  std::cout << "Best objective value: " << mcSol.get_weight() << std::endl;
  const std::vector<int>& solution = mcSol.get_assignments();
  for (int i=0; i < solution.size(); ++i) {
    std::cout << "Index " << i << " : " << solution[i] << std::endl;
  }
  return 0;
}
```

The constructor to `MaxCutInstance` takes a file name as an argument and loads the instance from the file. We construct `heur`, our heuristic, by passing it the instance, the runtime limit, an indicator for whether to run validity checks after the run is completed, and a pointer to an object specifying callbacks (which we pass as `NULL`). The heuristic is run in the constructor, so the run will be complete once this second line of code in `main` is done running.

The best solution found during a heuristic run can be accessed with the `get_best_solution` function, which returns a reference to a `MaxCutSimpleSolution` object. More detailed information about the history of best solution values found can be obtained via the `History` function, or by directly accessing the `past_solution_values_` and `past_solution_times_` vectors within the heuristic. A more detailed history of the solutions themselves can be accessed by setting `validation=true` when instantiating the heuristic, and then accessing the `past_solutions_` vector within the heuristic. The `get_weight` and `get_assignments` functions can be used to obtain the objective value and node assignments for this solution.

To build this code, we need to provide the `C++` compiler with the header files for the MQLib project, which can be done with the `-I` compiler flag. Further, we need to provide the library `bin/MQLib.a` to the compiler. Thus the compilation line would be `g++ -Iinclude linked.cpp bin/MQLib.a`, run from the main MQLib folder. You could indicate a different output executable name with `g++ -Iinclude -o linked.out linked.cpp bin/MQLib.a`.

Running the resulting executable with `./a.out` should yield output similar to:

```
Best objective value: 3
Index 0 : -1
Index 1 : 1
Index 2 : 1
```

### Running a Heuristic Given Its Name

A heuristic can be run based on its name using the `HeuristicFactory` class, which is used in the MQLib executable both to dispatch heuristics specified on the command line and also to dispatch heuristics in the hyper-heuristic.

The below code uses a `HeuristicFactory` object to run the QUBO heuristic `GLOVER1998a` on the QUBO problem instance [bin/sampleQUBO.txt](sampleQUBO.txt).

```
#include <iostream>
#include "heuristics/heuristic_factory.h"

int main(int argc, char** argv) {
  QUBOInstance qi("bin/sampleQUBO.txt");
  HeuristicFactory factory;
  QUBOHeuristic* heur = factory.RunQUBOHeuristic("GLOVER1998a", qi,
                                                 10.0, false, NULL);
  const QUBOSimpleSolution& sol = heur->get_best_solution();
  std::cout << "Best objective value: " << sol.get_weight() << std::endl;
  const std::vector<int>& solution = sol.get_assignments();
  for (int i=0; i < solution.size(); ++i) {
    std::cout << "Index " << i << " : " << solution[i] << std::endl;
  }
  delete heur;
  return 0;
}
```

Of note, the `HeuristicFactory` returns a pointer to a heuristic that was allocated with the `new` operator, and it is the responsibility of the calling code to free this memory once it is no longer needed using the `delete` operator (otherwise the code will have a memory leak).

The code can be compiled as with the first example, and provides similar output:

```
Best objective value: 6
Index 0 : 1
Index 1 : 0
Index 2 : 1
```

### Running a Heuristic for One Problem on an Instance of the Other

Due to the simple reductions available from Max-Cut to QUBO and from QUBO to Max-Cut, it is simple to obtain solutions for a problem instance for one type of problem using a heuristic for the other problem. For instance, consider using the hyper-heuristic provided in the code base (which is a `MaxCutHeuristic`) to solve a QUBO instance. This can be done by reducing the QUBO instance to a Max-Cut instance, solving using the hyper-heuristic, and then converting the resulting solution back to a solution for the original QUBO instance:

```
#include <iostream>
#include "heuristics/maxcut/hyperheuristic.h"

int main(int argc, char** argv) {
  QUBOInstance qi("bin/sampleQUBO.txt");
  MaxCutInstance mi(qi);
  std::string selectedHeuristic;
  MaxCutHyperheuristic heur(mi, 10.0, false, NULL, 144, &selectedHeuristic);
  std::cout << "Selected heuristic: " << selectedHeuristic << std::endl;

  const MaxCutSimpleSolution& mcSol = heur.get_best_solution();
  QUBOSimpleSolution quboSol(mcSol, qi, NULL);
  std::cout << "Best objective value: " << quboSol.get_weight() << std::endl;
  const std::vector<int>& solution = quboSol.get_assignments();
  for (int i=0; i < solution.size(); ++i) {
    std::cout << "Index " << i << " : " << solution[i] << std::endl;
  }
  return 0;
}
```

We construct the `MaxCutInstance` named `mi` out of the `QUBOInstance` named `qi` by passing `qi` as the sole argument to the `MaxCutInstance` constructor (a `QUBOInstance` can be obtained from a `MaxCutInstance` by flipping this construction). When constructing (and thus running) the hyper-heuristic, we pass two additional arguments that are not passed to other heuristics: a seed and a pointer to a string, which it populates with the name of the heuristic that was selected for the run. After the heuristic is done running, we convert the best solution from a `MaxCutSimpleSolution` named `mcSol` to a `QUBOSimpleSolution` for the original instance named `quboSol` by constructing a `QUBOSimpleSolution` using `mcSol`, the original instance `qi`, and a `NULL` pointer for the associated `QUBOHeuristic`.

Building the code as before, we obtain similar output, indicating the heuristic that was selected by the hyper-heuristic and the best solution to the original QUBO instance:

```
Selected heuristic: BEASLEY1998TS
Best objective value: 6
Index 0 : 1
Index 1 : 0
Index 2 : 1
```

### Programmatically Providing Problem Instances

In the examples thus far, new problem instances have been specified by loading an instance from a specified file. However, problem instances can also be constructed by passing a vector of edges (Max-Cut) or a vector of non-zero off-diagonal terms in the input matrix along with a vector containing the main diagonal of the input matrix (QUBO).

As an example, consider programmatically creating the QUBO matrix specified in the file [bin/sampleQUBO.txt](sampleQUBO.txt), which has main diagonal `[5, 3, 1]` and off-diagonal values `Q_12 = Q_21 = -6` and `Q_23 = Q_32 = -1`. We could do this and use it to run the `GLOVER1998a` heuristic with:

```
#include <iostream>
#include <vector>
#include "heuristics/qubo/glover1998a.h"

int main(int argc, char** argv) {
  const double m[] = {5.0, 3.0, 1.0};
  std::vector<double> mainDiag(m, m+3);
  std::vector<Instance::InstanceTuple> offDiag;
  offDiag.push_back(Instance::InstanceTuple(std::make_pair(1, 2), -6.0));
  offDiag.push_back(Instance::InstanceTuple(std::make_pair(2, 3), -1.0));
  QUBOInstance qi(offDiag, mainDiag, 3);
  Glover1998a heur(qi, 10.0, false, NULL);

  const QUBOSimpleSolution& sol = heur.get_best_solution();
  std::cout << "Best objective value: " << sol.get_weight() << std::endl;
  const std::vector<int>& solution = sol.get_assignments();
  for (int i=0; i < solution.size(); ++i) {
    std::cout << "Index " << i << " : " << solution[i] << std::endl;
  }
  return 0;
}
```

When constructing the off-diagonal entries of the instance, we create a vector of type `Instance::InstanceTuple`, which is a typedef of type `std::pair<std::pair<int, int>, double>`. Variable numbers are 1-indexed when providing this information. The output of this code is the same as the previous example:

```
Best objective value: 6
Index 0 : 1
Index 1 : 0
Index 2 : 1
```

A `MaxCutInstance` can be constructed similarly, by providing an edge list as a vector of type `Instance::InstanceTuple`. Node numbers are 1-indexed in the edge list. For instance, a Max-Cut instance with three nodes and a link of weight 3 between nodes 1 and 2 and a link of weight -2 between nodes 2 and 3 could be constructed with:

```
std::vector<Instance::InstanceTuple> edgeList;
edgeList.push_back(Instance::InstanceTuple(std::make_pair(1, 2), 3.0));
edgeList.push_back(Instance::InstanceTuple(std::make_pair(2, 3), -2.0));
MaxCutInstance mi(edgeList, 3);
```

### Accessing New Best Solutions and Implementing Custom Termination Criteria with Callbacks

All the examples thus far have provided a runtime limit that was used to determine when to stop running a heuristic. However, custom termination criteria are also available using the callback mechanism for heuristics. Heuristics routinely call the `Report` function to report their progress (regardless of whether they have a new best solution to report or not), and each call to `Report` generates a callback that provides the best solution encountered thus far, whether this best solution was just encountered (aka it has not been reported in previous callbacks), the elapsed runtime thus far, and the iteration count of the heuristic (if the heuristic reports this information). To receive callbacks, code running a Max-Cut heuristic should create a class that extends the `MaxCutCallback` class and code running a QUBO heuristic should create a class that extends the `QUBOCallback` class; a `Report` function of the callback class will be called whenever the heuristic calls `Report`. As an example, we can use the callback mechanism to run the `BURER2002` heuristic on the `bin/sampleMaxCut.txt` instance until it has not reported a new best solution for more than 5 seconds:

```
#include <iostream>
#include <vector>
#include "heuristics/maxcut/burer2002.h"

class Burer2002Callback : public MaxCutCallback {
 public:
  Burer2002Callback() :
    lastNewBest_(0.0) {}

  bool Report(const MaxCutSimpleSolution& sol, bool newBest,
               double runtime) {
    if (newBest) {
      std::cout << "New best " << sol.get_weight() << " at " << runtime <<
        std::endl;
      lastNewBest_ = runtime;
    }
    if (runtime - lastNewBest_ >= 5.0) {
      std::cout << "Exiting at runtime " << runtime << std::endl;
      return false;
    } else {
      return true;  // Keep going
    }
  }

  bool Report(const MaxCutSimpleSolution& sol, bool newBest,
              double runtime, int iter) {
    return Report(sol, newBest, runtime);
  }
    
 private:
  double lastNewBest_;  // Elapsed seconds when last new best sol. was found
};

int main(int argc, char** argv) {
  MaxCutInstance mi("bin/sampleMaxCut.txt");
  Burer2002Callback callback;
  Burer2002 heur(mi, 10.0, false, &callback);
  return 0;
}
```

Though a runtime limit is still provided in the constructor of the heuristic, this is ignored when running with a callback (provided as the last argument to the heuristic). Callbacks must handle both `Report` calls that include an iteration count and those that do not; in our case we simply called the non-iteration count version of the `Report` function when provided with an iteration count. This code should produce output similar to:

```
New best 3 at 6.4e-05
Exiting at runtime 5.00007
```
