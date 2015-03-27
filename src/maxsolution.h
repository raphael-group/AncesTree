/*
 *  maxsolution.h
 *
 *   Created on: 5-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef MAXSOLUTION_H
#define MAXSOLUTION_H

#include "utils.h"
#include "realmatrix.h"
#include "clonaltree.h"

namespace vaff {

class MaxSolution
{
public:
  struct Triple
  {
    RealMatrix _F;
    RealMatrix _U;
    ClonalTree _T;
    
    bool operator==(const Triple& other) const
    {
      return _U == other._U && _T == other._T && _F == other._F;
    }
    
    bool operator!=(const Triple& other) const
    {
      return !(operator==(other));
    }
    
    friend std::ostream& operator<<(std::ostream& out,
                                    const Triple& triple);
    
    friend std::istream& operator>>(std::istream& in,
                                    Triple& triple);
  };
  
  typedef std::vector<Triple> TripleVector;
  typedef TripleVector::const_iterator TripleVectorIt;
  
  MaxSolution(const RealMatrix& F);
  
  MaxSolution();
  
  void add(const Triple& triple)
  {
    _triples.push_back(triple);
  }
  
  const Triple& solution(int index) const
  {
    assert(0 <= index && index < _triples.size());
    return _triples[index];
  }
  
  int size() const
  {
    return _triples.size();
  }
  
  void clear()
  {
    _triples.clear();
  }
  
  bool present(const Triple& sol) const
  {
    for (TripleVectorIt it = _triples.begin(); it != _triples.end(); ++it)
    {
      if (*it == sol)
      {
        return true;
      }
    }
    
    return false;
  }
  
  const RealMatrix& getF() const
  {
    return _F;
  }
  
  void remapLabels(int max_cluster_size);
  
  double computeVafDelta(int sol_idx) const;
  
  void printVafDelta(int sol_idx, std::ostream& out) const;
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const MaxSolution& solution);
  
  friend std::istream& operator>>(std::istream& in,
                                  MaxSolution& solution);
  
private:
  RealMatrix _F;
  TripleVector _triples;
  
  typedef TripleVector::iterator TripleVectorNonConstIt;
};

} // namespace vaff

#endif // MAXSOLUTION_H