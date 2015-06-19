/*
 *  analyse.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "ancestrymatrix.h"
#include "readcountmatrix.h"
//#include "probancestrygraphyoshiko.h"
#include "probancestrygraph.h"
#include <cmath>

#include <fstream>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <READ_COUNTS> <ALPHA> <BETA> <GAMMA> where" << std::endl
      << "  <READ_COUNTS>      is the ancestry matrix file" << std::endl
      << "  <ALPHA>            alpha parameter (ancestry)" << std::endl
      << "  <BETA>             beta parameter (equality)" << std::endl
      << "  <GAMMA>            gamma parameter (CI)" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc != 5)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }
  
  double alpha = -1;
  sscanf(argv[2], "%lf", &alpha);
  if (!(0 <= alpha && alpha <= 0.5))
  {
    std::cerr << "Error: alpha must be in [0,0.5]" << std::endl;
    return 1;
  }
  
  double beta = -1;
  sscanf(argv[3], "%lf", &beta);
  if (!(0.5 <= beta && beta <= 1))
  {
    std::cerr << "Error: beta must be in [0.5,1]" << std::endl;
    return 1;
  }
  
  double gamma = -1;
  sscanf(argv[4], "%lf", &gamma);
  if (!(0 <= gamma && gamma <= 1))
  {
    std::cerr << "Error: gamma must be in [0,1]" << std::endl;
    return 1;
  }

  ReadCountMatrix R;
  std::string read_count_matrix = argv[1];
  
  if (read_count_matrix != "-")
  {
    std::ifstream in(read_count_matrix.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << read_count_matrix << "' for reading" << std::endl;
      return 1;
    }
    in >> R;
    in.close();
  }
  else
  {
    std::cin >> R;
  }
  
  AncestryMatrix A(R, 0);
  
  ProbAncestryGraph G(A, R, alpha, gamma);
  
  StlIntMatrix toOrginalColumns;
  G.removeCycles(A, alpha, toOrginalColumns);

  ReadCountMatrix newR = R.collapse(toOrginalColumns);
  newR.remapLabels(toOrginalColumns, R);
  RealIntervalMatrix CI;
  newR.computeConfidenceIntervals(CI, gamma);
  
  ProbAncestryGraph H;
  G.contract(A, toOrginalColumns, beta, H);

  std::cout << CI;
  std::cout << newR;
  
  return 0;
}