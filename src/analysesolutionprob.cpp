/*
 *  analysesolutionprob.cpp
 *
 *   Created on: 17-jan-2015
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "maxsolution.h"
#include "solutiongraph.h"
#include "ancestrymatrix.h"
#include "probancestrygraph.h"
#include "probcomparison.h"
#include <stdlib.h>
#include <fstream>
#include <iterator>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <SOLUTION_IDX> <SOLUTION> <ANCESTRY_MATRIX> <TYPE> where" << std::endl
      << "  <SOLUTION_IDX>     is the solution index" << std::endl
      << "  <SOLUTION>         is the solution file, specify '-' to use stdin" << std::endl
      << "  <ANCESTRY_MATRIX>  is the ancestry matrix file" << std::endl
      << "  <TYPE>             0 : delta VAF\n"
         "                     1 : clustered\n"
         "                     2 : ancestral\n"
         "                     3 : incomparable\n"
         "                     4 : clustered (median)\n"
         "                     5 : ancestral (median)\n"
         "                     6 : incomparable (median)\n"
         "                     7 : delta VAF (mean)" << std::endl;
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
  
  if (!(0 <= sol_idx && sol_idx < solution.size()))
  {
    std::cerr << "Invalid solution index " << sol_idx << "; it must be in the range [0, " << solution.size() << ")" << std::endl;
    return 1;
  }
  
  AncestryMatrix A;
  std::string filename_ancestry_matrix = argv[3];
  
  if (filename_ancestry_matrix != "-")
  {
    std::ifstream in(filename_ancestry_matrix.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << filename_ancestry_matrix << "' for reading" << std::endl;
      return 1;
    }
    in >> A;
    in.close();
  }
  
  const MaxSolution::Triple& sol = solution.solution(sol_idx);
  int type = atoi(argv[4]);
  switch(type)
  {
    case 0:
      solution.printVafDelta(sol_idx, std::cout);
      break;
    case 1:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.clustered();
        std::copy(S.begin(), S.end(),
                  std::ostream_iterator<double>(std::cout, "\n"));
      }
      break;
    case 2:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.ancestral();
        std::copy(S.begin(), S.end(),
                  std::ostream_iterator<double>(std::cout, "\n"));
      }
      break;
    case 3:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.incomparable();
        std::copy(S.begin(), S.end(),
                  std::ostream_iterator<double>(std::cout, "\n"));
      }
      break;
    case 4:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.clustered();
        if (S.size() > 0)
          std::cout << comp.median(S);
        std::cout << std::endl;
      }
      break;
    case 5:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.ancestral();
        if (S.size() > 0)
          std::cout << comp.median(S);
        std::cout << std::endl;
      }
      break;
    case 6:
      {
        StlDoubleVector S;
        ProbComparison comp(A, sol);
        S = comp.incomparable();
        if (S.size() > 0)
          std::cout << comp.median(S);
        std::cout << std::endl;
      }
      break;
    case 7:
      std::cout << solution.computeVafDelta(sol_idx) << std::endl;
      break;
    default:
      std::cerr << "Invalid type" << std::endl;
      return 1;
  }
  
  return 0;
}
