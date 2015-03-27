/*
 *  probcomparison.h
 *
 *   Created on: 26-mar-2015
 *       Author: M. El-Kebir
 */

#ifndef PROBCOMPARISON_H
#define PROBCOMPARISON_H

#include "utils.h"
#include "maxsolution.h"
#include "clonaltree.h"

#include <set>
#include <algorithm>

namespace vaff {
  
class ProbComparison
{
public:
  ProbComparison(const AncestryMatrix& A,
                 const MaxSolution::Triple& solution);
  
  StlDoubleVector clustered() const;
  StlDoubleVector ancestral() const;
  StlDoubleVector incomparable() const;
  
  static double median(StlDoubleVector& S)
  {
    if (S.size() % 2 == 1)
    {
      std::nth_element(S.begin(), S.begin() + S.size() / 2, S.end());
      return S[S.size() / 2];
    }
    else
    {
      std::nth_element(S.begin(), S.begin() + S.size() / 2 - 1, S.end());
      std::nth_element(S.begin(), S.begin() + S.size() / 2, S.end());
      
      return (S[S.size() / 2] + S[S.size() / 2 - 1]) / 2;
    }
  }
  
private:
  typedef std::pair<int, int> StlIntPair;
  typedef std::set<StlIntPair> StlIntPairSet;
  
  const AncestryMatrix& _A;
  const MaxSolution::Triple& _solution;
  StlIntMatrix _toMutationsFromSol;
  StlIntVector _toSolCluster;
  StlIntPairSet _clustered;
  StlIntPairSet _ancestral;
  StlIntPairSet _incomparable;
  
  void constructMappings();
  
  void determinePairs();
};

} // namespace vaff

#endif // PROBCOMPARISON_H
