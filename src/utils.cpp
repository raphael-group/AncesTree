/*
 *  utils.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include <sstream>
#include <string>

namespace vaff {
  
lemon::Tolerance<double> g_tol(1e-6);
  
std::ostream& operator<<(std::ostream& out, const StlBoolMatrix& M)
{
  int m = M.size();
  int n = M.empty() ? -1 : M[0].size();
  
  out << m << std::endl;
  out << n << std::endl;
  
  for (int i = 0; i < m; ++i)
  {
    for (int k = 0; k < n; ++k)
    {
      out << M[i][k] << " ";
    }
    out << std::endl;
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, StlBoolMatrix& M)
{
  int m = -1, n = -1;
  
  std::string line;
  std::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  std::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  M = StlBoolMatrix(m, StlBoolVector(n, false));
  for (int i = 0; i < m; ++i)
  {
    std::getline(in, line);
    ss.clear();
    ss.str(line);

    for (int j = 0; j < n; ++j)
    {
      int b = -1;
      ss >> b;
      
      if (!(b == 0 || b == 1))
      {
        throw std::runtime_error("Invalid entry");
      }
      
      M[i][j] = b;
    }
  }
  
  return in;
}
  
std::ostream& operator<<(std::ostream& out, const StlDoubleMatrix& M)
{
  int m = M.size();
  int n = M.empty() ? -1 : M[0].size();
  
  out << m << std::endl;
  out << n << std::endl;
  
  for (int i = 0; i < m; ++i)
  {
    for (int k = 0; k < n; ++k)
    {
      out << M[i][k] << " ";
    }
    out << std::endl;
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, StlDoubleMatrix& M)
{
  int m = -1, n = -1;
  
  std::string line;
  std::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  std::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  M = StlDoubleMatrix(m, StlDoubleVector(n, 0));
  for (int i = 0; i < m; ++i)
  {
    std::getline(in, line);
    ss.clear();
    ss.str(line);
    
    for (int j = 0; j < n; ++j)
    {
      ss >> M[i][j];
    }
  }
  
  return in;
}
  
std::ostream& operator<<(std::ostream& out, const StlRealIntervalMatrix& M)
{
  int m = M.size();
  int n = M.empty() ? -1 : M[0].size();
  
  out << m << std::endl;
  out << n << std::endl;
  
  for (int i = 0; i < m; ++i)
  {
    for (int k = 0; k < n; ++k)
    {
      out << M[i][k].first << " ";
    }
    out << std::endl;
  }
  
  std::cout << std::endl;
  
  out << m << std::endl;
  out << n << std::endl;
  
  for (int i = 0; i < m; ++i)
  {
    for (int k = 0; k < n; ++k)
    {
      out << M[i][k].second << " ";
    }
    out << std::endl;
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, StlRealIntervalMatrix& M)
{
  int m = -1, n = -1;
  
  std::string line;
  std::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  std::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  M = StlRealIntervalMatrix(m, StlRealIntervalVector(n));
  for (int i = 0; i < m; ++i)
  {
    std::getline(in, line);
    ss.clear();
    ss.str(line);
    
    for (int j = 0; j < n; ++j)
    {
      ss >> M[i][j].first;
    }
  }
  
  // skip blank line
  std::getline(in, line);
  
  int m2 = -1, n2 = -1;
  
  std::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> m2;
  
  if (m2 != m)
  {
    throw std::runtime_error("Error: m and m' should match");
  }
  
  std::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n2;
  
  if (n2 != n)
  {
    throw std::runtime_error("Error: n and n' should match");
  }
  
  for (int i = 0; i < m; ++i)
  {
    std::getline(in, line);
    ss.clear();
    ss.str(line);
    
    for (int j = 0; j < n; ++j)
    {
      ss >> M[i][j].second;
    }
  }
  
  return in;
}
  
StlBoolVector discretize(const StlDoubleVector& v)
{
  int n = v.size();
  StlBoolVector res(n, false);
  for (int i = 0; i < n; ++i)
  {
    if (v[i] != 0)
    {
      res[i] = true;
    }
  }
  return res;
}
  
StlBoolMatrix discretize(const StlDoubleMatrix& M)
{
  int m = M.size();
  int n = m != 0 ? M.front().size() : 0;

  StlBoolMatrix res(m, StlBoolVector(n, false));
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      if (M[i][j] != 0)
      {
        res[i][j] = true;
      }
    }
  }
  
  return res;
}
  
} // namespace vaff
