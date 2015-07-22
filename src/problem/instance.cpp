#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "problem/instance.h"

void Instance::AddLink(int n1, int n2, double weight,
		       std::vector<std::vector<std::pair<int, double> > >* links,
		       std::vector<InstanceTuple>* all,
		       std::vector<double>* selfLinks, bool selfLinkAsError) {
  if (n1 == n2) {
    if (selfLinkAsError) {
      std::cout << "Self-link encountered" << std::endl;
      exit(1);
    }
    if (selfLinks) {  // Only store self-link data if selfLinks is non-null
      (*selfLinks)[n1] = weight;
    }
    return;
  }

  (*links)[n1].push_back(std::pair<int, double>(n2, weight));
  (*links)[n2].push_back(std::pair<int, double>(n1, weight));
  if (n2 > n1) {
    all->push_back(InstanceTuple(std::make_pair(n1, n2), weight));
  } else {
    all->push_back(InstanceTuple(std::make_pair(n2, n1), weight));
  }
}

void Instance::Load(int dimension,
                    const std::vector<InstanceTuple>& provided,
                    std::vector<std::vector<std::pair<int, double> > >* links,
                    std::vector<InstanceTuple>* all,
                    std::vector<double>* selfLinks, bool selfLinkAsError) {
  if (!links || !all) {
    std::cout << "Invalid pointers passed to Instance::Load" << std::endl;
    exit(1);
  }
  links->clear();
  all->clear();
  if (selfLinks) {
    selfLinks->clear();
  }
  if (dimension <= 0) {
    std::cout << "Illegal dimension: " << dimension << std::endl;
    exit(1);
  }

  // Initialize state based on stated dimension
  for (int ct=0; ct < dimension; ++ct) {
    links->push_back(std::vector<std::pair<int, double> >());
    if (selfLinks) {
      selfLinks->push_back(0.0);
    }
  }

  for (auto iter=provided.begin(); iter != provided.end(); ++iter) {
    if (iter->first.first < 1 || iter->first.first > dimension) {
      std::cout << "Illegal first node in tuple (nodes are 1-indexed): " <<
        iter->first.first << std::endl;
      exit(1);
    }
    if (iter->first.second < 1 || iter->first.second > dimension) {
      std::cout << "Illegal second node in tuple (nodes are 1-indexed): " <<
        iter->first.second << std::endl;
      exit(1);
    }
    AddLink(iter->first.first-1, iter->first.second-1, iter->second, links, all,
            selfLinks, selfLinkAsError);
  }
}

void Instance::Load(const std::string& filename,
		    std::vector<std::vector<std::pair<int, double> > >* links,
		    std::vector<InstanceTuple>* all,
		    std::vector<double>* selfLinks, bool selfLinkAsError) {
  if (!links || !all) {
    std::cout << "Invalid pointers passed to Instance::Load" << std::endl;
    exit(1);
  }
  links->clear();
  all->clear();
  if (selfLinks) {
    selfLinks->clear();
  }
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    std::cout << "File cannot be opened: " << filename << std::endl;
    exit(1);    
  }
  std::string line;
  int dimension;
  int numLines;
  bool readDimension = false;
  int dataLinesRead = 0;
  while (getline(file, line)) {
    if (line.length() == 0 || line.at(0) == '#') {
      continue;  // Ignore empty lines or comments
    }
    if (!readDimension) {
      // First line in file contains "dimension" and "numLines", space-separated
      if (sscanf(line.c_str(), "%d %d", &dimension, &numLines) != 2) {
        std::cout << "Illegal first line: " << line << std::endl;
        exit(1);
      }
      if (dimension <= 0) {
        std::cout << "Illegal dimension: " << dimension << std::endl;
        exit(1);
      }
      if (numLines < 0) {
        std::cout << "Illegal number of data lines: " << numLines << std::endl;
        exit(1);
      }

      // Initialize state based on stated dimension
      for (int ct=0; ct < dimension; ++ct) {
	links->push_back(std::vector<std::pair<int, double> >());
	if (selfLinks) {
	  selfLinks->push_back(0.0);
	}
      }
      readDimension = true;
    } else {
      if (dataLinesRead == numLines) {
        std::cout << "Extra data line: " << line << std::endl;
        exit(1);
      }
      int n1, n2;
      double weight;
      if (sscanf(line.c_str(), "%d %d %lf", &n1, &n2, &weight) != 3) {
        std::cout << "Illegal data line: " << line << std::endl;
        exit(1);
      }
      if (n1 < 1 || n1 > dimension) {
        std::cout << "Illegal first node in data line (nodes are 1-indexed): "
                  << line << std::endl;
        exit(1);
      }
      if (n2 < 1 || n2 > dimension) {
        std::cout << "Illegal second node in data line (nodes are 1-indexed): "
                  << line << std::endl;
        exit(1);
      }
      AddLink(n1-1, n2-1, weight, links, all, selfLinks, selfLinkAsError);
      ++dataLinesRead;
    }
  }
  if (file.bad()) {
    std::cout << "IO error reading file " << filename << std::endl;
    exit(1);
  }
  if (dataLinesRead < numLines) {
    std::cout << "Not enough data lines in " << filename << std::endl;
    exit(1);
  }
}
