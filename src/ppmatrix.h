/*
 *  ppmatrix.h
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#ifndef PPMATRIX_H
#define PPMATRIX_H

#include "matrix.h"
#include "utils.h"

namespace vaff {
  
class PpMatrix : public Matrix
{
public:
  PpMatrix();
  
  PpMatrix(int m, int n);
  
  const StlBoolMatrix& getMatrix() const
  {
    return _C;
  }
  
  bool operator()(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _C[row][col];
  }
  
  void set(int row, int col, bool val)
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    _C[row][col] = val;
  }
  
  int rowSum(int row) const
  {
    assert(0 <= row && row < _m);
    
    int res = 0;
    for (int j = 0; j < _n; ++j)
    {
      if (_C[row][j])
      {
        ++res;
      }
    }
    
    return res;
  }
  
  int isRowContained(int row1, int row2) const
  {
    assert(0 <= row1 && row1 < _m);
    assert(0 <= row2 && row2 < _m);
    
    int diff = 0;
    for (int j = 0; j < _n; ++j)
    {
      int val1 = _C[row1][j] ? 1 : 0;
      int val2 = _C[row2][j] ? 1 : 0;
      
      if (val1 == 1 && val2 == 0)
      {
        return -1;
      }
      else if (val1 == 0 && val2 == 1)
      {
        ++diff;
      }
    }
    
    return diff;
  }
  
  bool isConflictFree() const;
  
  bool hasConflict(int col1, int col2) const;
  
  bool isContained(int col1, int col2) const;
  
  bool isDisjoint(int col1, int col2) const;
  
  bool operator==(const PpMatrix& other) const
  {
    return _m == other._m && _n == other._n && _C == other._C;
  }
  
  bool operator!=(const PpMatrix& other) const
  {
    return !this->operator==(other);
  }
  
  unsigned long twoPowerColSum(int col) const
  {
    assert(0 <= col && col < _n);
    
    unsigned long res = 0;
    for (int i = 0; i < _m; ++i)
    {
      if (_C[i][col])
      {
        res |= (1 << i);
      }
    }
    
    return res;
  }
  
  double prob(int j) const
  {
    assert(0 <= j && j < _m);
    return _prob[j];
  }
  
  void setProb(int j, double p)
  {
    assert(0 <= j && j < _m);
    _prob[j] = p;
  }
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const PpMatrix& matrix);
  
  friend std::istream& operator>>(std::istream& in,
                                  PpMatrix& matrix);
  
protected:
  StlBoolMatrix _C;
  StlDoubleVector _prob;
};
  
} // namespace vaff

#endif // PPMATRIX_H