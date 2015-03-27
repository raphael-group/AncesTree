/*
 *  ancestrymatrix.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "ancestrymatrix.h"
#include <cmath>

namespace vaff {

AncestryMatrix::AncestryMatrix()
  : _n(0)
  , _C()
{
}
  
AncestryMatrix::AncestryMatrix(const ReadCountMatrix& R,
                               int order)
  : _n(R.getNrRows())
  , _C(_n, StlDoubleVector(_n, 0))
{
  StlDoubleVector log_fact;
  
  int max_count = 0;
  const int m = R.getNrCols();
  const int n = R.getNrRows();
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    max_count = std::max(max_count, std::max(R.getAlt(j, i), R.getRef(j, i)));
  }
  
  constructLogFactorialTable(4*max_count, log_fact);
  
  for (int p = 0; p < n; ++p)
  {
    for (int q = 0; q < n; ++q)
    {
      _C[p][q] = prob(R, order, log_fact, p, q);
    }
  }
}

void AncestryMatrix::constructLogFactorialTable(const int n,
                                                StlDoubleVector& log_fact)
{
  log_fact = StlDoubleVector(n+1, 0);
  log_fact[0] = 0;
  for (int i = 1; i <= n; ++i)
  {
    log_fact[i] = log_fact[i - 1] + log(i);
  }
}
  
double AncestryMatrix::prob(const ReadCountMatrix& R,
                            int order,
                            const StlDoubleVector& log_fact,
                            int p,
                            int q)
{
  assert(0 <= p && p < R.getNrRows());
  assert(0 <= q && q < R.getNrRows());
  
  const int m = R.getNrCols();
  StlDoubleVector prob_vector(m);
  for (int i = 0; i < m; ++i)
  {
//    std::cerr << "(" << i << "," << R.getRowLabel(p) << "," << R.getRowLabel(q) << ") : " << prob(R, log_fact, i, p, q) << std::endl;
//    std::cerr << "(" << i << "," << p << "," << q << ") : " << prob(R, log_fact, i, p, q) << std::endl;

    prob_vector[i] = std::max(0.0, std::min(1.0, prob(R, log_fact, i, p, q)));
  }
  std::sort(prob_vector.begin(), prob_vector.end());
  return prob_vector[order];
}
  
double AncestryMatrix::prob(const ReadCountMatrix& R,
                            const StlDoubleVector& log_fact,
                            int i,
                            int p,
                            int q)
{
  assert(0 <= p && p < R.getNrRows());
  assert(0 <= q && q < R.getNrRows());
  assert(0 <= i && i < R.getNrCols());
  
  int alt_p = R.getAlt(p, i);
  int ref_p = R.getRef(p, i);
  int alt_q = R.getAlt(q, i);
  int ref_q = R.getRef(q, i);
  
  // we want to compute the probability that p is an ancestor to q
  if (alt_p == 0 && ref_p == 0 && alt_q == 0 && ref_q == 0)
  {
    return 1;
  }
  else if (alt_p == 0 && ref_p == 0)
  {
    return 0;
  }
  else if (alt_q == 0 && ref_q == 0)
  {
    return 1;
  }
  else
  {
    return 1 - g(log_fact, alt_q+1, ref_q+1, alt_p+1, ref_p+1);
  }
}
  
double AncestryMatrix::g(const StlDoubleVector& log_fact,
                         int a, int b, int c, int d)
{
  int aa = std::min(a, c);
  int bb = std::min(b, d);
  int cc = aa;
  int dd = bb;

  double res = 0.5;
  
  while (aa < a)
  {
    res += h(log_fact, aa, bb, cc, dd) / aa;
    ++aa;
  }

  while (bb < b)
  {
    res -= h(log_fact, aa, bb, cc, dd) / bb;
    ++bb;
  }
  
  while (cc < c)
  {
    res -= h(log_fact, aa, bb, cc, dd) / cc;
    ++cc;
  }
  
  while (dd < d)
  {
    res += h(log_fact, aa, bb, cc, dd) / dd;
    ++dd;
  }
  
  assert(res != INFINITY);
  assert(res != -INFINITY);
  assert(res != NAN);
  
  return res;
}
  
double AncestryMatrix::h(const StlDoubleVector& log_fact,
                         int a, int b, int c, int d)
{
  return exp(log_beta(log_fact, a+c, b+d) - (log_beta(log_fact, a, b) + log_beta(log_fact, c, d)));
}
  
double AncestryMatrix::log_beta(const StlDoubleVector& log_fact,
                                int x, int y)
{
  assert(x + y - 1 < log_fact.size());
  return log_fact[x - 1] + log_fact[y - 1] - log_fact[x + y - 1];
}
  
std::ostream& operator<<(std::ostream& out,
                         const AncestryMatrix& matrix)
{
  out << matrix._C;
  return out;
}

std::istream& operator>>(std::istream& in,
                         AncestryMatrix& matrix)
{
  in >> matrix._C;
  matrix._n = matrix._C.size();
  return in;
}
  
}; // namespace vaff