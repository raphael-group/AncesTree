/*
 *  realintervalmatrix.cpp
 *
 *   Created on: 7-jan-2015
 *       Author: M. El-Kebir
 */

#include "realintervalmatrix.h"

namespace vaff {

RealIntervalMatrix::RealIntervalMatrix()
  : Matrix()
  , _C()
{
}
  
RealIntervalMatrix::RealIntervalMatrix(int m, int n)
  : Matrix(m, n)
  , _C(m, StlRealIntervalVector(n))
{
}
  
void RealIntervalMatrix::obtainPointEstimates(RealMatrix& F) const
{
  F = RealMatrix(_m, _n);
  for (int i = 0; i < _m; ++i)
  {
    for (int j = 0; j < _n; ++j)
    {
      const RealInterval& interval = _C[i][j];
      F.set(i, j, 0.5 * (interval.first + interval.second));
    }
  }
}
  
std::ostream& operator<<(std::ostream& out,
                         const RealIntervalMatrix& matrix)
{
  out << matrix._C;
  if (matrix._labelsSet)
  {
    out << std::endl;
    
    for (int i = 0; i < matrix._m; ++i)
    {
      out << matrix._rowLabel[i] << " ";
    }
    out << std::endl;
    
    for (int j = 0; j < matrix._n; ++j)
    {
      out << matrix._colLabel[j] << " ";
    }
    out << std::endl;
  }
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         RealIntervalMatrix& matrix)
{
  in >> matrix._C;
  
  matrix._m = matrix._C.size();
  matrix._n = matrix._C.empty() ? 0 : matrix._C.front().size();
  matrix._rowLabel.resize(matrix._m);
  matrix._colLabel.resize(matrix._n);
  
  return in;
}
    
} // namespace vaff
