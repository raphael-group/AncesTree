/*
 *  realintervalmatrix.h
 *
 *   Created on: 7-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef REALINTERVALMATRIX_H
#define REALINTERVALMATRIX_H

#include "matrix.h"
#include "utils.h"
#include "realmatrix.h"

namespace vaff {
  
class RealIntervalMatrix : public Matrix
{
public:
  RealIntervalMatrix();
  
  RealIntervalMatrix(int m, int n);
  
  const StlRealIntervalMatrix& getMatrix() const
  {
    return _C;
  }
  
  const RealInterval& operator()(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _C[row][col];
  }
  
  void set(int row, int col, RealInterval interval)
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    _C[row][col] = interval;
  }
   
  bool operator==(const RealIntervalMatrix& other) const
  {
    return _m == other._m && _n == other._n && _C == other._C;
  }
  
  bool operator!=(const RealIntervalMatrix& other) const
  {
    return !this->operator==(other);
  }
  
  double meanIntervalLength(bool excludeTight) const
  {
    double mean = 0;
    int count = 0;
    
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        if (!excludeTight || _C[i][j].first != _C[i][j].second)
        {
          mean += _C[i][j].second - _C[i][j].first;
          ++count;
        }
      }
    }
    
    return mean / count;
  }
  
  int tightAndZeroIntervalCount() const
  {
    int count = 0;
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        if (_C[i][j].first == 0 && _C[i][j].second == 0)
        {
          ++count;
        }
      }
    }
    return count;
  }
  
  int tightIntervalCount() const
  {
    int count = 0;
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        if (_C[i][j].first == _C[i][j].second)
        {
          ++count;
        }
      }
    }
    return count;
  }
  
  void obtainPointEstimates(RealMatrix& F) const;
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const RealIntervalMatrix& matrix);
  
  friend std::istream& operator>>(std::istream& in,
                                  RealIntervalMatrix& matrix);
  
protected:
  StlRealIntervalMatrix _C;
};
  
}

#endif // REALINTERVALMATRIX_H