/*
 *  readcountmatrix.h
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef READCOUNTMATRIX_H
#define READCOUNTMATRIX_H

#include "matrix.h"
#include "utils.h"
#include "realintervalmatrix.h"

namespace vaff {
  
class ReadCountMatrix : public Matrix
{
public:
  ReadCountMatrix();
  
  ReadCountMatrix(int m, int n);
  
  const StlIntMatrix& getC() const
  {
    // alternate
    return _C;
  }
  
  const StlIntMatrix& getD() const
  {
    // reference
    return _D;
  }
  
  int getAlt(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _C[row][col];
  }
  
  int getRef(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _D[row][col];
  }
  
  IntPair operator()(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return std::make_pair(_C[row][col], _D[row][col]);
  }
  
  void set(int row, int col, int alt, int ref)
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    _C[row][col] = alt;
    _D[row][col] = ref;
  }
  
  double coverage() const
  {
    double res = 0;
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        res += _C[i][j] + _D[i][j];
      }
    }
    return res / (_m * _n);
  }
  
  void remapLabels(const StlIntMatrix& toOrgCols,
                   const ReadCountMatrix& orgR,
                   int max_cluster_size = 10);
  
  void computeConfidenceIntervals(RealIntervalMatrix& CI, double gamma) const;
  
  void computePointEstimates(RealMatrix& F) const;
  
  ReadCountMatrix collapse(StlIntMatrix& toOrgColumns) const;
  
  bool operator==(const ReadCountMatrix& other) const
  {
    return _m == other._m && _n == other._n && _C == other._C && _D == other._D;
  }
  
  bool operator!=(const ReadCountMatrix& other) const
  {
    return !this->operator==(other);
  }
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const ReadCountMatrix& matrix);
  
  friend std::istream& operator>>(std::istream& in,
                                  ReadCountMatrix& matrix);
  
protected:
  StlIntMatrix _C;
  StlIntMatrix _D;
};

} // namespace vaff

#endif // READCOUNTMATRIX_H
