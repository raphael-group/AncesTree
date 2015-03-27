/*
 *  ancestrymatrix.h
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef ANCESTRYMATRIX_H
#define ANCESTRYMATRIX_H

#include "realmatrix.h"
#include "readcountmatrix.h"
#include "utils.h"
#include <lemon/tolerance.h>

namespace vaff {
  
class AncestryMatrix
{
public:
  AncestryMatrix();
  
  AncestryMatrix(const ReadCountMatrix& R,
                 int order);
  
  double operator()(int row, int col) const
  {
    assert(0 <= row && row < _n);
    assert(0 <= col && col < _n);
    
    return _C[row][col];
  }
  
  int getNrRows() const
  {
    return _n;
  }
  
  int getNrCols() const
  {
    return _n;
  }
  
  void writeAntiSymmetricElements(const ReadCountMatrix& R,
                                  const StlIntMatrix& M,
                                  std::ostream& out) const
  {
    for (StlIntMatrixIt it = M.begin(); it != M.end(); ++it)
    {
      writeAntiSymmetricElements(R, *it, out);
    }
  }
  
  void writeAntiSymmetricElements(const ReadCountMatrix& R,
                                  const StlIntVector& S,
                                  std::ostream& out) const
  {
    StlDoubleVector log_fact;
    constructLogFactorialTable(100000, log_fact);
    
    int m = R.getNrCols();
    lemon::Tolerance<double> tol(1e-3);
    lemon::Tolerance<double> tol2(1e-2);
    for (StlIntVectorIt it1 = S.begin(); it1 != S.end(); ++it1)
    {
      int i = *it1;
      for (StlIntVectorIt it2 = it1; it2 != S.end(); ++it2)
      {
        int j = *it2;
        if (i == j)
        {
          continue;
        }
//        if (!tol.nonZero(1 - (_C[i][j] + _C[j][i])))
        {
          out << R.getRowLabel(i);
          for (int k = 0; k < m; ++k)
          {
            out << "\t" << R.getAlt(i, k) << "\t" << R.getRef(i, k);
          }
          out << std::endl;
          
          out << R.getRowLabel(j);
          for (int k = 0; k < m; ++k)
          {
            out << "\t" << R.getAlt(j, k) << "\t" << R.getRef(j, k);
          }
          out << std::endl;
          
          out << "\t";
          int count = 0;
          double min_i_j = 1;
          for (int k = 0; k < m; ++k)
          {
            double p = std::max(0., AncestryMatrix::prob(R, log_fact, k, i, j));
            min_i_j = std::min(min_i_j, p);
            if (tol2.nonZero(p) && tol2.less(p, 1))
            {
              ++count;
            }
            out << "\t\t" << p;
          }
          out << std::endl;
          out << "\t";
          double min_j_i = 1;
          for (int k = 0; k < m; ++k)
          {
            double p = std::max(0., AncestryMatrix::prob(R, log_fact, k, j, i));
            min_j_i = std::min(min_j_i, p);
            out << "\t\t" << p;
          }
          out << std::endl;
          
          if (std::max(min_i_j, min_j_i) < 0.05)
          {
            std::cout << "-";
          }
          else if (std::max(min_i_j, min_j_i) > 0.9)
          {
            std::cout << "+";
          }
          
          if (!tol.nonZero(1 - (_C[i][j] + _C[j][i])))
          {
            std::cout << "*";
          }
          
          std::cout << count << std::endl;
          
          out << std::endl;
        }
      }
    }
  }
  
  void writeAntiSymmetricElements(const ReadCountMatrix& R,
                                  std::ostream& out) const
  {
    StlDoubleVector log_fact;
    constructLogFactorialTable(100000, log_fact);
    
    int m = R.getNrCols();
    lemon::Tolerance<double> tol(1e-3);
    lemon::Tolerance<double> tol2(1e-2);
    for (int i = 0; i < _n; ++i)
    {
      for (int j = i + 1; j < _n; ++j)
      {
//        if (!tol.nonZero(1 - (_C[i][j] + _C[j][i])))
        {
          out << R.getRowLabel(i);
          for (int k = 0; k < m; ++k)
          {
            out << "\t" << R.getAlt(i, k) << "\t" << R.getRef(i, k);
          }
          out << std::endl;
          
          out << R.getRowLabel(j);
          for (int k = 0; k < m; ++k)
          {
            out << "\t" << R.getAlt(j, k) << "\t" << R.getRef(j, k);
          }
          out << std::endl;

          out << "\t";
          int count = 0;
          for (int k = 0; k < m; ++k)
          {
            double p = std::max(0., AncestryMatrix::prob(R, log_fact, k, i, j));
            if (tol2.nonZero(p) && tol2.less(p, 1))
            {
              ++count;
            }
            out << "\t\t" << p;
          }
          out << std::endl;
          out << "\t";
          for (int k = 0; k < m; ++k)
          {
            double p = std::max(0., AncestryMatrix::prob(R, log_fact, k, j, i));
            out << "\t\t" << p;
          }
          out << std::endl;
          if (count > 1)
            std::cout << "*" << std::endl;

          out << std::endl;
        }
      }
    }
  }
  
  void antiSymmetricElements(const StlIntMatrix& M,
                             int& antiPairs,
                             int& totalPairs) const
  {
    antiPairs = totalPairs = 0;
    for (StlIntMatrixIt it = M.begin(); it != M.end(); ++it)
    {
      const StlIntVector& S = *it;
      totalPairs += S.size() * S.size();
      antiPairs += antiSymmetricElements(S);
    }
  }
  
  int antiSymmetricElements(const StlIntVector& S) const
  {
    int res = 0;
    lemon::Tolerance<double> tol(1e-3);
    
    for (StlIntVectorIt it1 = S.begin(); it1 != S.end(); ++it1)
    {
      int j = *it1;
      for (StlIntVectorIt it2 = it1; it2 != S.end(); ++it2)
      {
        int k = *it2;
        if (j == k)
          continue;
        
        if (!tol.nonZero(1 - (_C[j][k] + _C[k][j])))
        {
          ++res;
        }
      }
    }
    
    return res*2 + S.size();
  }
  
  int antiSymmetricElements() const
  {
    int res = 0;
    lemon::Tolerance<double> tol(1e-3);
    for (int i = 0; i < _n; ++i)
    {
      for (int j = i + 1; j < _n; ++j)
      {
        if (!tol.nonZero(1 - (_C[i][j] + _C[j][i])))
        {
          ++res;
        }
      }
    }
    return res*2 + _n;
  }
  
  static void constructLogFactorialTable(const int n,
                                         StlDoubleVector& log_fact);
  
  static double prob(const ReadCountMatrix& R,
                     int order,
                     const StlDoubleVector& log_fact,
                     int p, int q);
  
  static double prob(const ReadCountMatrix& R,
                     const StlDoubleVector& log_fact,
                     int i, int p, int q);
  
  static double g(const StlDoubleVector& log_fact,
                  int a, int b, int c, int d);
  
  static double h(const StlDoubleVector& log_fact,
                  int a, int b, int c, int d);
  
  static double log_beta(const StlDoubleVector& log_fact,
                         int x, int y);
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const AncestryMatrix& matrix);
  
  friend std::istream& operator>>(std::istream& in,
                                  AncestryMatrix& matrix);

private:
  int _n;
  StlDoubleMatrix _C;
};
  
} // namespace vaff

#endif // ANCESTRYMATRIX_H
