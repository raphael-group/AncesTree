/*
 *  analysesolution.cpp
 *
 *   Created on: 11-mar-2015
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "maxsolution.h"
#include "solutiongraph.h"
#include "ancestrymatrix.h"
#include "probancestrygraph.h"
#include "comparison.h"
#include <stdlib.h>
#include <fstream>
#include <cstring>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <SOLUTION> <REFERENCE_SOLUTION> <WHITELIST_SOLUTION> where" << std::endl
      << "  <SOLUTION>            is the solution file (use '-' for stdin)" << std::endl
      << "  <REFERENCE_SOLUTION>  is the reference solution" << std::endl
      << "  <WHITELIST_SOLUTION>  is the solution to be used as a whitelist" << std::endl;
  out << "Output:\n"
         "  coverage\n"
         "  delta F\n"
         "  delta U\n"
         "  recall clustered\n"
         "  recall ancestral\n"
         "  recall incomparable\n"
         "  accuracy clustered\n"
         "  accuracy ancestral\n"
         "  accuracy incomparable" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }
  
  MaxSolution solution;
  {
    if (strcmp(argv[1], "-"))
    {
      std::ifstream in(argv[1]);
      if (!in.good())
      {
        std::cerr << "Error: failed to open '" << argv[1] << "' for reading" << std::endl;
        return 1;
      }
      in >> solution;
      in.close();
    }
    else
    {
      std::cin >> solution;
    }
  }
  
  MaxSolution refSolution;
  {
    std::ifstream in(argv[2]);
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << argv[2] << "' for reading" << std::endl;
      return 1;
    }
    in >> refSolution;
    in.close();
  }
  
  MaxSolution whitelistSolution;
  {
    std::ifstream in(argv[3]);
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << argv[3] << "' for reading" << std::endl;
      return 1;
    }
    in >> whitelistSolution;
    in.close();
  }
  
  Comparison comparison(refSolution, solution.solution(0), whitelistSolution.solution(0));
  
  double coverage = comparison.coverage();
  
  double deltaF = comparison.deltaF();
  double deltaU = comparison.deltaU();
  
  double clustered, ancestral, incomparable;
  comparison.recallB(clustered, ancestral, incomparable);
  
  double clustered_a, ancestral_a, incomparable_a;
  comparison.accuracyB(clustered_a, ancestral_a, incomparable_a);
  
  std::cout << coverage << "\t"
            << deltaF << "\t"
            << deltaU << "\t"
            << clustered << "\t"
            << ancestral << "\t"
            << incomparable << "\t"
            << clustered_a << "\t"
            << ancestral_a << "\t"
            << incomparable_a << std::endl;
  
  return 0;
}