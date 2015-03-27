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
  out << "Usage: " << argv0 << " <READ_COUNTS> <ANCESTRY_MATRIX> <ALPHA> <BETA> <GAMMA> <TYPE> (<IDX>) where" << std::endl
      << "  <READ_COUNTS>      is the ancestry matrix file" << std::endl
      << "  <ANCESTRY_MATRIX>  is the read count matrix file" << std::endl
      << "  <ALPHA>            alpha parameter (ancestry)" << std::endl
      << "  <BETA>             beta parameter (equality)" << std::endl
      << "  <GAMMA>            gamma parameter (CI)" << std::endl
      << "  <TYPE>             0 : #nodes in-deg 0, #nodes CI > 0.5, largest arborescence size\n"
         "                     1 : intra cluster coherence\n"
         "                     2 : inter cluster coherence\n"
         "                     3 : is G a DAG?\n"
         "                     4 : #non-trivial SCC\n"
         "                     5 : is G's transitive closure G itself\n"
         "                     6 : is H a DAG?\n"
         "                     7 : sizes of SCCs of G\n"
         "                     8 : ancestry subgraph of G of specified SCC\n"
         "                     9 : fraction of (clustered) mutation pairs summing to 1\n"

         "                     10: fraction of mutations pairs of specified SCC summing to 1\n"
         "                     11: read data for mutations summing to 1\n"
         "                     12: read data for clustered mutations summing to 1\n"
         "                     13: read data for clustered mutations of specified SCC summing to 1\n"
         "                     14: write yoshiko input\n"
         "                     15: ancestry graph G\n"
         "                     16: coverage\n"
         "                     17: ancestry graph H\n"
         "                     18: sample alpha and output #clusters\n"
         "                     19: sample beta (given alpha) and output #edges" << std::endl
      << "  <IDX>              SCC index" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc != 7 && argc != 8)
  {
    printUsage(argv[0], std::cerr);
    return 1;
  }
  
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
  
  AncestryMatrix A;
  std::string filename_ancestry_matrix = argv[2];
  
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
  
  //ProbAncestryGraphYoshiko G(A, R, alpha, gamma);
  ProbAncestryGraph G(A, R, alpha, gamma);
  
  StlIntMatrix toOrginalColumns;
  G.removeCycles(A, alpha, toOrginalColumns);

  ReadCountMatrix newR = R.collapse(toOrginalColumns);
  newR.remapLabels(toOrginalColumns, R);
  RealIntervalMatrix CI;
  newR.computeConfidenceIntervals(CI, gamma);
  
  ProbAncestryGraph H;
  G.contract(A, toOrginalColumns, beta, H);
//  H.writeDOT(newR, CI, std::cout);
  
  int type = atoi(argv[6]);
  switch(type)
  {
    case 0:
      std::cout << H.numberOfNodesInDeg0() << "\t"
                << G.numberOfNodesInfCI() << "\t"
                << H.largestArborescence() << std::endl;
      break;
    case 1:
      {
        ProbAncestryGraph::StringVector label;
        ProbAncestryGraph::createLabels(toOrginalColumns, label);
        
        ProbAncestryGraph::printIntraClusterCoherence(H.getG(),
                                                      H.getNodeToColumnMap(),
                                                      label,
                                                      A, std::cout);
        
        ProbAncestryGraph::printClusterSize(toOrginalColumns, std::cerr);
      }
      break;
    case 2:
      {
        ProbAncestryGraph::StringVector label;
        ProbAncestryGraph::createLabels(toOrginalColumns, label);
        
        ProbAncestryGraph::printInterClusterCoherence(H.getG(),
                                                      H.getNodeToColumnMap(),
                                                      label,
                                                      A, std::cout);
      }
      break;
    case 3:
//      std::cout << G.numberOfNodesInfCI() << "\t" << G.isDAG() << std::endl;
      std::cout << G.isDAG() << std::endl;
      break;
    case 4:
      std::cout << G.numberOfNonTrivialSCC() << std::endl;
      break;
    case 5:
      std::cout << G.isTransitive() << std::endl;
      break;
    case 6:
      std::cout << H.isDAG() << std::endl;
      break;
    case 7:
      {
        int i = 0;
        for (StlIntMatrixIt it = toOrginalColumns.begin(); it != toOrginalColumns.end(); ++it, ++i)
        {
          std::cout << i << ": " << it->size() << std::endl;
        }
      }
      break;
    case 8:
      if (argc != 8)
      {
        printUsage(argv[0], std::cerr);
        return 1;
      }
      else
      {
        int idx = atoi(argv[7]);
        if (0 <= idx && idx < toOrginalColumns.size())
        {
          G.writeDOT(R, toOrginalColumns[idx], std::cout);
        }
        else
        {
          std::cerr << "Invalid SCC index" << std::endl;
          return 1;
        }
      }
      break;
    case 9:
      {
        int n_square = A.getNrRows() * A.getNrRows();
        int anti_count = A.antiSymmetricElements();
        double f = (float)anti_count / (float)n_square;
        std::cout << f << ",";
      }
      {
        int totalCount;
        int antiCount;
        A.antiSymmetricElements(toOrginalColumns, antiCount, totalCount);
        double f = (float)antiCount / (float)totalCount;
        std::cout << f << std::endl;
      }
      break;
    case 10:
      if (argc != 8)
      {
        printUsage(argv[0], std::cerr);
        return 1;
      }
      else
      {
        int idx = atoi(argv[7]);
        if (0 <= idx && idx < toOrginalColumns.size())
        {
          int n_square = toOrginalColumns[idx].size();
          n_square *= n_square;
          int anti_count = A.antiSymmetricElements(toOrginalColumns[idx]);
          double f = (float)anti_count / (float)n_square;
          std::cout << anti_count << "/" << n_square << " = " << f << std::endl;
        }
        else
        {
          std::cerr << "Invalid SCC index" << std::endl;
          return 1;
        }
      }
      break;
    case 11:
      A.writeAntiSymmetricElements(R, std::cout);
      break;
    case 12:
      A.writeAntiSymmetricElements(R, toOrginalColumns, std::cout);
    case 13:
      if (argc != 8)
      {
        printUsage(argv[0], std::cerr);
        return 1;
      }
      else
      {
        int idx = atoi(argv[7]);
        if (0 <= idx && idx < toOrginalColumns.size())
        {
          A.writeAntiSymmetricElements(R, toOrginalColumns[idx], std::cout);
        }
        else
        {
          std::cerr << "Invalid SCC index" << std::endl;
          return 1;
        }
      }
      break;
    case 14:
//      G.writeYoshiko(A, R, alpha, std::cout);
      break;
    case 15:
      G.writeDOT(R, toOrginalColumns, std::cout);
      break;
    case 16:
      std::cout << R.coverage() << std::endl;
      break;
    case 17:
      H.writeDOT(newR, std::cout);
      H.isDAG();
      break;
    case 18:
      for (double a = 0; a <= 0.5; a += 0.01)
      {
        ProbAncestryGraph GG(A, R, a, gamma);
        
        StlIntMatrix toOrginalColumns2;
        GG.removeCycles(A, a, toOrginalColumns2);
        
        std::cout << "," << toOrginalColumns2.size();
      }
      std::cout << std::endl;
      break;
    case 19:
      for (double b = 0.5; b <= 1; b += 0.01)
      {
        ProbAncestryGraph HH;
        G.contract(A, toOrginalColumns, b, HH);
        
        std::cout << "," << lemon::countArcs(HH.getG());
      }
      std::cout << std::endl;
      break;
    default:
      std::cerr << "Invalid type" << std::endl;
      return 1;
  }
  
  return 0;
}