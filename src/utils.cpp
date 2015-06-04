/*
 *  utils.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include <sstream>
#include <string>
#include <cstdio>       // EOF

namespace vaff {
  
lemon::Tolerance<double> g_tol(1e-6);
  
inline std::istream& getline(std::istream& is, std::string& t)
{
  // copied from: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
  t.clear();
  
  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  
  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();
  
  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if(sb->sgetc() == '\n')
          sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if(t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}
  
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
  vaff::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  vaff::getline(in, line);
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
    vaff::getline(in, line);
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
  vaff::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  vaff::getline(in, line);
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
    vaff::getline(in, line);
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
  vaff::getline(in, line);
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  vaff::getline(in, line);
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
    vaff::getline(in, line);
    ss.clear();
    ss.str(line);
    
    for (int j = 0; j < n; ++j)
    {
      ss >> M[i][j].first;
    }
  }
  
  // skip blank line
  vaff::getline(in, line);
  
  int m2 = -1, n2 = -1;
  
  vaff::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> m2;
  
  if (m2 != m)
  {
    throw std::runtime_error("Error: m and m' should match");
  }
  
  vaff::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n2;
  
  if (n2 != n)
  {
    throw std::runtime_error("Error: n and n' should match");
  }
  
  for (int i = 0; i < m; ++i)
  {
    vaff::getline(in, line);
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
