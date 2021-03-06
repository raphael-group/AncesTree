/*
 *  machinainput.cpp
 *
 *   Created on: 30-nov-2017
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "maxsolution.h"
#include "solutiongraph.h"
#include <stdlib.h>
#include <fstream>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <SOLUTION_IDX> <SOLUTION> <THRESHOLD> <CLUSTER_SIZE> <BETA> where" << std::endl
      << "  <SOLUTION_IDX>  is the solution index" << std::endl
      << "  <SOLUTION>      is the solution file, specify '-' to use stdin" << std::endl
      << "  <THRESHOLD>     usage threshold"  << std::endl
      << "  <BETA>          beta" << std::endl;
}


int main(int argc, char** argv)
{
  if (argc != 5)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }
  
  int sol_idx = atoi(argv[1]);
  std::string filename = argv[2];
  
  double threshold = -1;
  sscanf(argv[3], "%lf", &threshold);
  if (!(0 <= threshold && threshold <= 1.0))
  {
    std::cerr << "Error: threshold must be in [0,1.0]" << std::endl;
    return 1;
  }
  
  double beta = -1;
  sscanf(argv[4], "%lf", &beta);
  if (!(0.5 <= beta && beta <= 1))
  {
    std::cerr << "Error: beta must be in [0.5,1]" << std::endl;
    return 1;
  }
  
  MaxSolution solution;
  
  if (filename != "-")
  {
    std::ifstream in(filename.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << argv[2] << "' for reading" << std::endl;
      return 1;
    }
    in >> solution;
    in.close();
  }
  else
  {
    std::cin >> solution;
  }
  
  solution.remapLabels(-1);
  SolutionGraph graph(solution.solution(sol_idx), threshold, beta);
  graph.writeEdgeList(std::cout);
  graph.writeLeaves(std::cerr);
  
  return 0;
}
