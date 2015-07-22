#ifndef HEURISTICS_HEURISTIC_FACTORY_H_
#define HEURISTICS_HEURISTIC_FACTORY_H_

#include <map>
#include <string>
class MaxCutHeuristic;
#include "problem/max_cut_heuristic.h"
#include "problem/max_cut_instance.h"
class QUBOHeuristic;
#include "problem/qubo_heuristic.h"
#include "problem/qubo_instance.h"

class HeuristicFactory {
 public:
  // Constructor, which creates function pointers stored in max_cut_map_ and
  // qubo_map that can be used to generate each heuristic.
  HeuristicFactory();


  // Run a heuristic (runs are performed in the heuristics' constructors). After
  // the run is complete, the function will return with a pointer to that
  // heuristic. If heuristic code is invalid, returns NULL.
  // NOTE: To avoid memory leaks, you are responsible for calling delete on the
  // pointer after you are done using it.
  MaxCutHeuristic* RunMaxCutHeuristic(const std::string& code,
                                      const MaxCutInstance& mi,
                                      double runtime_limit,
                                      bool validation,
                                      MaxCutCallback *mc);
  QUBOHeuristic* RunQUBOHeuristic(const std::string& code,
                                  const QUBOInstance& qi,
                                  double runtime_limit,
                                  bool validation,
                                  QUBOCallback *qc);
  
  // Checks if the passed heuristic code is valid.
  bool ValidMaxCutHeuristicCode(const std::string& code) const {
    return max_cut_map_.find(code) != max_cut_map_.end();
  }
  bool ValidQUBOHeuristicCode(const std::string& code) const {
    return qubo_map_.find(code) != qubo_map_.end();
  }

  // Get all heuristic codes, in alphabetical order
  void MaxCutHeuristicCodes(std::vector<std::string>* codes);
  void QUBOHeuristicCodes(std::vector<std::string>* codes);

  // Get the instance (it will be converted from a provided instance if it was
  // not provided and has not already been generated via conversion). Recipient
  // is not responsible for freeing memory from the result (and should not call
  // delete on the received object). If factory was constructed without providing
  // any instances, function will return NULL.
  MaxCutInstance* GetMaxCutInstance();
  QUBOInstance* GetQUBOInstance();

  // Output the codes and descriptions for each heuristic.
  void PrintHeuristicCodes();

 private:  
  // Mapping from heuristic names to creator functions
  typedef std::pair<MaxCutHeuristic*(*)(const MaxCutInstance& mi,
                                        double runtime_limit,
                                        bool validation, MaxCutCallback *mc),
    std::string> MaxCutCreator;
  std::map<std::string, MaxCutCreator> max_cut_map_;
  typedef std::pair<QUBOHeuristic*(*)(const QUBOInstance& mi,
                                      double runtime_limit,
                                      bool validation, QUBOCallback *qc),
    std::string> QUBOCreator;
  std::map<std::string, QUBOCreator> qubo_map_;
};

#endif
