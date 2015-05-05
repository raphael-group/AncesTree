/*
 *  ppmatrix.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "ppmatrix.h"
#include <sstream>

namespace vaff {

PpMatrix::PpMatrix()
  : Matrix()
  , _C()
  , _prob()
{
}
  
PpMatrix::PpMatrix(int m, int n)
  : Matrix(m, n)
  , _C(m, StlBoolVector(n, false))
  , _prob(m, 1)
{
}
  
bool PpMatrix::isConflictFree() const
{
  for (int col1 = 0; col1 < _n; ++col1)
  {
    for (int col2 = col1 + 1; col2 < _n; ++col2)
    {
      if (hasConflict(col1, col2))
      {
        return false;
      }
    }
  }
  
  return true;
}

bool PpMatrix::hasConflict(int col1, int col2) const
{
  assert(0 <= col1 && col1 < _n);
  assert(0 <= col2 && col2 < _n);
  bool one_one = false, one_zero = false, zero_one = false;

  for (int i = 0; i < _m; ++i)
  {
    one_one |= _C[i][col1] && _C[i][col2];
    one_zero |= _C[i][col1] && !_C[i][col2];
    zero_one |= !_C[i][col1] && _C[i][col2];
  }
  
  return one_one && one_zero && zero_one;
}

bool PpMatrix::isContained(int col1, int col2) const
{
  assert(0 <= col1 && col1 < _n);
  assert(0 <= col2 && col2 < _n);

  for (int i = 0; i < _m; ++i)
  {
    if (_C[i][col1] && !_C[i][col2])
      return false;
  }
  
  return true;
}

bool PpMatrix::isDisjoint(int col1, int col2) const
{
  assert(0 <= col1 && col1 < _n);
  assert(0 <= col2 && col2 < _n);
  
  for (int i = 0; i < _m; ++i)
  {
    if (_C[i][col1] && _C[i][col2])
      return false;
  }
  
  return true;
}

std::ostream& operator<<(std::ostream& out,
                         const PpMatrix& matrix)
{
  out << matrix._C;
  out << std::endl;
  
  for (int j = 0; j < matrix._m; ++j)
  {
    out << matrix._prob[j] << " ";
  }
  out << std::endl;
  
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         PpMatrix& matrix)
{
  in >> matrix._C;
  
  matrix._m = matrix._C.size();
  matrix._n = matrix._C.empty() ? 0 : matrix._C.front().size();
  matrix._rowLabel.resize(matrix._m);
  matrix._colLabel.resize(matrix._n);
  
  matrix._prob = StlDoubleVector(matrix._m, 1);
  
  std::string line;
  vaff::getline(in, line);
  vaff::getline(in, line);
  std::stringstream ss(line);
  for (int j = 0; j < matrix._m; ++j)
  {
    ss >> matrix._prob[j];
  }
  
  return in;
}

} // namespace vaff
