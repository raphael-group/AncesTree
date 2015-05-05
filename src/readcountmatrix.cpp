/*
 *  readcountmatrix.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "readcountmatrix.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/beta.hpp>

namespace vaff {

ReadCountMatrix::ReadCountMatrix()
  : Matrix()
  , _C()
  , _D()
{
}
  
ReadCountMatrix::ReadCountMatrix(int m, int n)
  : Matrix(m, n)
  , _C(m, StlIntVector(n, 0))
  , _D(m, StlIntVector(n, 0))
{
}
  
void ReadCountMatrix::remapLabels(const StlIntMatrix& toOrgColumns,
                                  const ReadCountMatrix& orgR,
                                  int max_cluster_size)
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  char buf[1024];
  for (int j = 0; j < _m; ++j)
  {
    const StlIntVector& orgColumns = toOrgColumns[j];
    std::string new_label = "";
    bool first = true;
    int count = 0;
    for (StlIntVectorIt it = orgColumns.begin(); it != orgColumns.end() && count < max_cluster_size; ++it, ++count)
    {
      int org_j = *it;
      if (first)
      {
        first = false;
      }
      else
      {
        new_label += "\n";
      }
      StringVector s2 ;
      boost::split(s2, orgR.getRowLabel(org_j), boost::is_any_of(","));
      new_label += s2.front();
    }
    if (count == max_cluster_size && orgColumns.size() > max_cluster_size)
    {
      snprintf(buf, 1024, "\n[%d more]", (int)orgColumns.size() - max_cluster_size);
      new_label += buf;
    }
    setRowLabel(j, new_label);
  }
}
  
std::ostream& operator<<(std::ostream& out,
                         const ReadCountMatrix& matrix)
{
  out << "gene_id";
  for (int j = 0; j < matrix._n; ++j)
  {
    out << "\t" << matrix.getColLabel(j) << "\t" << matrix.getColLabel(j);
  }
  out << std::endl;
  
  for (int i = 0; i < matrix._m; ++i)
  {
    out << matrix.getRowLabel(i);
    for (int j = 0; j < matrix._n; ++j)
    {
      out << "\t" << matrix._D[i][j] << "\t" << matrix._C[i][j];
    }
    out << std::endl;
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         ReadCountMatrix& matrix)
{
  typedef std::vector<std::string> StringVector;

  std::string line;
  std::getline(in, line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t"));
  
  if (s.empty())
  {
    throw std::runtime_error("Error: empty sample labels");
  }
  
  // pop first element
  s = StringVector(++s.begin(), s.end());
  
  if (s.size() % 2 != 0)
  {
    throw std::runtime_error("Error: odd number of samples");
  }
  
  matrix._n = s.size() / 2;
  matrix._colLabel = StringVector(matrix._n);
  for (int j = 0; j < matrix._n; ++j)
  {
    if (s[2*j] != s[2*j + 1])
    {
      throw std::runtime_error("Error: unequal sample label between ref and alt");
    }
    matrix._colLabel[j] = s[2*j];
  }
  
  matrix._rowLabel.clear();
  matrix._C.clear();
  matrix._D.clear();
  while (std::getline(in, line).good() && line != "")
  {
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != 1+ 2*matrix._n)
    {
      throw std::runtime_error("Error: invalid number of columns");
    }
    
    matrix._rowLabel.push_back(s[0]);
    matrix._C.push_back(StlIntVector(matrix._n));
    matrix._D.push_back(StlIntVector(matrix._n));
    
    // pop first element
    s = StringVector(++s.begin(), s.end());
    
    for (int j = 0; j < matrix._n; ++j)
    {
      int ref = -1;
      int alt = -1;
      
      if (s[2*j] == "")
      {
        ref = 0;
      }
      else
      {
        ref = atoi(s[2*j].c_str());
        if (ref < 0)
        {
          throw std::runtime_error("Error: ref count is negative");
        }
      }
      if (s[2*j+1] == "")
      {
        alt = 0;
      }
      else
      {
        alt = atoi(s[2*j+1].c_str());
        if (ref < 0)
        {
          throw std::runtime_error("Error: alt count is negative");
        }
      }
      
      matrix._C.back()[j] = alt;
      matrix._D.back()[j] = ref;
    }
  }
  
  matrix._m = matrix._rowLabel.size();
  
  return in;
}
  
ReadCountMatrix ReadCountMatrix::collapse(StlIntMatrix& toOrgColumns) const
{
  int nrNewMutations = toOrgColumns.size();
  ReadCountMatrix R(nrNewMutations, _n);
  R._colLabel = _colLabel;
  char buf[1024];
  
  for (int i = 0; i < nrNewMutations; ++i)
  {
    snprintf(buf, 1024, "%d", i);
    R._rowLabel[i] = buf;
    const StlIntVector& M = toOrgColumns[i];
    for (StlIntVectorIt it = M.begin(); it != M.end(); ++it)
    {
      for (int sample = 0; sample < _n; ++sample)
      {
        R._C[i][sample] += _C[*it][sample];
        R._D[i][sample] += _D[*it][sample];
      }
    }
  }
  
  return R;
}
  
void ReadCountMatrix::computeConfidenceIntervals(RealIntervalMatrix& CI,
                                                 double gamma) const
{
  typedef boost::math::beta_distribution<double> BetaDistribution;
  
  CI = RealIntervalMatrix(_n, _m);
  
  // samples
  for (int i = 0; i < _n; ++i)
  {
    // mutations
    for (int j = 0; j < _m; ++j)
    {
      if (_C[j][i] == 0 && _D[j][i] == 0)
      {
        // degenerate case:
        CI.set(i, j, RealInterval(0, 0));
      }
      else
      {
        BetaDistribution B(1 + getAlt(j, i), 1 + getRef(j, i));
        double UB;
        if (_D[j][i] == 0)
        {
          UB = 1;
        }
        else
        {
          UB = boost::math::quantile(B, 1 - gamma / 2);
        }
        
        double LB;
        if (_C[j][i] == 0)
        {
          LB = 0;
        }
        else
        {
          LB = boost::math::quantile(B, gamma / 2);
        }
        
        assert(0 <= LB && LB <= UB && UB <= 1);
        CI.set(i, j, RealInterval(LB, UB));
      }
    }
  }
}
  
void ReadCountMatrix::computePointEstimates(RealMatrix& F) const
{
  F = RealMatrix(_n, _m);
  for (int i = 0; i < _n; ++i)
  {
    F.setRowLabel(i, getColLabel(i));
    for (int j = 0; j < _m; ++j)
    {
      double ref = getRef(j, i);
      double alt = getAlt(j, i);
      if (ref + alt != 0)
        F.set(i, j,  alt / (ref + alt));
    }
  }
  
  for (int j = 0; j < _m; ++j)
  {
    F.setColLabel(j, getRowLabel(j));
  }
}
  
} // namespace vaff
