#ifndef PROBLEM_INSTANCE_H_
#define PROBLEM_INSTANCE_H_

#include <string>
#include <utility>
#include <vector>

class Instance {
 public:
  typedef std::pair<std::pair<int, int>, double> InstanceTuple;

  static void Load(int dimension,
                   const std::vector<InstanceTuple>& provided,
                   std::vector<std::vector<std::pair<int, double> > >* links,
		   std::vector<InstanceTuple>* all,
                   std::vector<double>* selfLinks, bool selfLinkAsError);

  static void Load(const std::string& filename,
		   std::vector<std::vector<std::pair<int, double> > >* links,
		   std::vector<InstanceTuple>* all,
		   std::vector<double>* selfLinks, bool selfLinkAsError);

 private:
  static void AddLink(int n1, int n2, double weight,
		      std::vector<std::vector<std::pair<int, double> > >* links,
		      std::vector<InstanceTuple>* all,
		      std::vector<double>* selfLinks, bool selfLinkAsError);
};

#endif
