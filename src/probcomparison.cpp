/*
 *  probcomparison.cpp
 *
 *   Created on: 26-mar-2015
 *       Author: M. El-Kebir
 */

#include "probcomparison.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace vaff {

ProbComparison::ProbComparison(const AncestryMatrix& A,
                               const MaxSolution::Triple& solution)
  : _A(A)
  , _solution(solution)
  , _toMutationsFromSol()
  , _toSolCluster()
{
  constructMappings();
  determinePairs();
}
  
void ProbComparison::constructMappings()
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  int n = _solution._F.getNrCols();
  
  _toMutationsFromSol = StlIntMatrix(n, StlIntVector());
  _toSolCluster = StlIntVector(_A.getNrRows(), -1);
  
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, _solution._F.getColLabel(j), boost::is_any_of(";"));
    
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      if (*it2 != "")
      {
        int org_j = boost::lexical_cast<int>(*it2);
        
        _toMutationsFromSol[j].push_back(org_j);
        _toSolCluster[org_j] = j;
      }
    }
  }
}
  
void ProbComparison::determinePairs()
{
  int org_j = 0;
  for (StlIntVectorIt it_1 = _toSolCluster.begin();
       it_1 != _toSolCluster.end(); ++it_1, ++org_j)
  {
    int j = *it_1;
    if (j == -1) continue;
    
    int org_k = org_j + 1;
    for (StlIntVectorIt it_2 = it_1 + 1;
         it_2 != _toSolCluster.end(); ++it_2, ++org_k)
    {
      int k = *it_2;
      if (k == -1) continue;
      
      if (j == k)
      {
        _clustered.insert(std::make_pair(std::min(org_j, org_k),
                                         std::max(org_j, org_k)));
      }
      else
      {
        if (_solution._T.ancestral(j, k))
        {
          _ancestral.insert(std::make_pair(org_j, org_k));
        }
        else if (_solution._T.ancestral(k, j))
        {
          _ancestral.insert(std::make_pair(org_k, org_j));
        }
        else
        {
          assert(_solution._T.incomparable(j, k));
          _incomparable.insert(std::make_pair(std::min(org_j, org_k),
                                              std::max(org_j, org_k)));
        }
      }
    }
  }
}

StlDoubleVector ProbComparison::clustered() const
{
  typedef StlIntPairSet::const_iterator StlIntPairSetIt;
  
  StlDoubleVector res;
  for (StlIntPairSetIt it = _clustered.begin();
       it != _clustered.end(); ++it)
  {
    res.push_back(_A(it->first, it->second));
    res.push_back(_A(it->second, it->first));
  }
  return res;
}

StlDoubleVector ProbComparison::ancestral() const
{
  typedef StlIntPairSet::const_iterator StlIntPairSetIt;
  
  StlDoubleVector res;
  for (StlIntPairSetIt it = _ancestral.begin();
       it != _ancestral.end(); ++it)
  {
    res.push_back(_A(it->first, it->second));
  }
  return res;
}
  
StlDoubleVector ProbComparison::incomparable() const
{
  typedef StlIntPairSet::const_iterator StlIntPairSetIt;
  
  StlDoubleVector res;
  for (StlIntPairSetIt it = _incomparable.begin();
       it != _incomparable.end(); ++it)
  {
    res.push_back(_A(it->first, it->second));
    res.push_back(_A(it->second, it->first));
  }
  return res;
}

} // namespace vaff