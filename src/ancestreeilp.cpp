/*
 *  ancestreeilp.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "utils.h"

#include "realmatrix.h"
#include "realintervalmatrix.h"
#include "maxsolution.h"
#include "intmaxilpsolver.h"
#include "clonaltree.h"
#include "ancestrymatrix.h"
//#include "probancestrygraphyoshiko.h"
#include "probancestrygraph.h"

#include <fstream>

using namespace vaff;

void printUsage(const char* argv0, std::ostream& out)
{
  out << "Usage: " << argv0 << " <READ_COUNTS> <ANCESTRY_MATRIX> <ALPHA> <BETA> <TIMELIMIT> where" << std::endl
      << "  <READ_COUNTS>      unclustered point estimates" << std::endl
      << "  <ANCESTRY_MATRIX>  clustered confidence intervals" << std::endl
      << "  <ALPHA>            alpha parameter (equality)" << std::endl
      << "  <BETA>             beta parameter (ancestry)" << std::endl
      << "  <GAMMA>            gamma parameter (CI)" << std::endl
      << "  <TIMELIMIT>        time limit in seconds (use -1 to disable time limit)" << std::endl;
}

bool readCountMatrix(const std::string& filename,
                     ReadCountMatrix& R)
{
  if (filename != "-")
  {
    std::ifstream in(filename.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << filename << "' for reading" << std::endl;
      return false;
    }
    in >> R;
    in.close();
  }
  
  return true;
}

bool readAncestryMatrix(const std::string& filename,
                        AncestryMatrix& A)
{
  if (filename != "-")
  {
    std::ifstream in(filename.c_str());
    if (!in.good())
    {
      std::cerr << "Error: failed to open '" << filename << "' for reading" << std::endl;
      return false;
    }
    in >> A;
    in.close();
  }
  
  return true;
}

bool readUnclustered(const std::string& filename,
                     RealMatrix& F_point_unclustered)
{
  std::ifstream in(filename.c_str());
  if (!in.good())
  {
    std::cerr << "Error: failed to open '" << filename << "' for reading" << std::endl;
    return false;
  }
  
  in >> F_point_unclustered;
  F_point_unclustered.setLabels(in);
  in.close();
  
  return true;
}

bool readClustered(const std::string& filename,
                   const int unclustered_n,
                   RealIntervalMatrix& F_interval_clustered,
                   StlIntMatrix& toUnclusteredColumns)
{
  std::ifstream in(filename.c_str());
  if (!in.good())
  {
    std::cerr << "Error: failed to open '" << filename << "' for reading" << std::endl;
    return false;
  }
  
  in >> F_interval_clustered;
  std::string line;
  vaff::getline(in, line);
  
  const int n = F_interval_clustered.getNrCols();
  toUnclusteredColumns = StlIntMatrix(n);
  for (int i = 0; i < n; ++i)
  {
    vaff::getline(in, line);
    std::stringstream ss(line);
    while (ss.good())
    {
      int val = -1;
      ss >> val;
      if (val == -1)
      {
        break;
      }
      if (!(0 <= val && val < unclustered_n))
      {
        std::cerr << "Error: column index " << val << " is out of bounds" << std::endl;
      }
      toUnclusteredColumns[i].push_back(val);
    }
  }
  
  in.close();
  
  return true;
}

int main(int argc, char** argv)
{
  if (argc != 7)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }

  int timeLimit = -1;
  timeLimit = atoi(argv[6]);
  
  double alpha = -1;
  sscanf(argv[3], "%lf", &alpha);
  if (!(0 <= alpha && alpha <= 0.5))
  {
    std::cerr << "Error: alpha must be in [0,0.5]" << std::endl;
    return 1;
  }
  
  double beta = -1;
  sscanf(argv[4], "%lf", &beta);
  if (!(0.5 <= beta && beta <= 1))
  {
    std::cerr << "Error: beta must be in [0.5,1]" << std::endl;
    return 1;
  }
  
  double gamma = -1;
  sscanf(argv[5], "%lf", &gamma);
  if (!(0 <= gamma && gamma <= 1))
  {
    std::cerr << "Error: gamma must be in [0,1]" << std::endl;
    return 1;
  }
  
  ReadCountMatrix R;
  if (!readCountMatrix(argv[1], R))
  {
    return 1;
  }
  
  AncestryMatrix A;
  if (!readAncestryMatrix(argv[2], A))
  {
    return 1;
  }
  
//  ProbAncestryGraphYoshiko G(A, R, alpha, gamma);
  ProbAncestryGraph G(A, R, alpha, gamma);
  
  StlIntMatrix toOrginalColumns;
  G.removeCycles(A, alpha, toOrginalColumns);
  
  RealMatrix F;
  R.computePointEstimates(F);
  
  ReadCountMatrix newR = R.collapse(toOrginalColumns);
  
  RealIntervalMatrix CI(newR.getNrCols(), newR.getNrRows());
  newR.computeConfidenceIntervals(CI, gamma);
  
  ProbAncestryGraph H;
  G.contract(A, toOrginalColumns, beta, H);
  
  std::cerr << "|V| = " << lemon::countNodes(H.getG()) << std::endl;
  std::cerr << "|A| = " << lemon::countArcs(H.getG()) << std::endl;
//  std::cerr << H.isDAG() << std::endl;
  
  IntMaxIlpSolver solver(H,
                         CI,
                         F,
                         toOrginalColumns,
                         timeLimit);
  
//  solver.exportModel("model.lp");
  MaxSolution solution(F);
  solver.solve(solution);
  
  std::cout << solution;
  
  return 0;
}
