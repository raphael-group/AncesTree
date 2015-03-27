/*
 *  print_max_solution.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "readcountmatrix.h"
#include "ancestrymatrix.h"
#include <stdlib.h>
#include <fstream>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <READ_COUNT_MATRIX> <ORDER> where" << std::endl
      << "  <READ_COUNT_MATRIX>  is the input file containing read counts\n"
      << "  <ORDER>              0 for minimum" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }
  
  std::string filename = argv[1];
  
  ReadCountMatrix R;
  if (filename != "-")
  {
    std::ifstream in(filename.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << argv[1] << "' for reading" << std::endl;
      return 1;
    }
    in >> R;
    in.close();
  }
  else
  {
    std::cin >> R;
  }
  
  int order = atoi(argv[2]);
  if (!(0 <= order && order < R.getNrCols()))
  {
    std::cerr << "Error: order has to be in [0, " << R.getNrCols() << "]" << std::endl;
    return 1;
  }
  
  AncestryMatrix M(R, order);
  std::cout << M;
  
  return 0;
}