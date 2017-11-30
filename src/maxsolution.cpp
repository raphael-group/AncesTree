/*
 *  maxsolution.cpp
 *
 *   Created on: 5-jan-2015
 *       Author: M. El-Kebir
 */

#include "maxsolution.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace vaff {

std::ostream& operator<<(std::ostream& out,
                         const MaxSolution::Triple& triple)
{
  out << triple._U;
  out << triple._T;
  out << triple._F;
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         MaxSolution::Triple& triple)
{
  in >> triple._U;
  in >> triple._T;
  in >> triple._F;
  triple._F.setLabels(in);
  
  return in;
}

MaxSolution::MaxSolution(const RealMatrix& F)
  : _F(F)
  , _triples()
{
}
  
MaxSolution::MaxSolution()
  : _F()
  , _triples()
{
}
  
double MaxSolution::computeVafDelta(int sol_idx) const
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  const Triple& sol = _triples[sol_idx];
  
  double delta = 0;
    
  int n = sol._F.getNrCols();
  int m = sol._F.getNrRows();
  int entries = 0;
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, sol._F.getColLabel(j), boost::is_any_of(";"));
    
    std::string new_label = "";
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      if (*it2 != "")
      {
        int org_j = boost::lexical_cast<int>(*it2);
        
        for (int i = 0; i < m; ++i)
        {
          delta += fabs(sol._F(i, j) - _F(i, org_j));
          ++entries;
        }
      }
    }
  }
  
  return delta / entries;
}
  
void MaxSolution::printVafDelta(int sol_idx, std::ostream& out) const
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  const Triple& sol = _triples[sol_idx];
  
  int n = sol._F.getNrCols();
  int m = sol._F.getNrRows();
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, sol._F.getColLabel(j), boost::is_any_of(";"));
    
    std::string new_label = "";
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      if (*it2 != "")
      {
        int org_j = boost::lexical_cast<int>(*it2);
        
        for (int i = 0; i < m; ++i)
        {
          out << fabs(sol._F(i, j) - _F(i, org_j)) << std::endl;
        }
      }
    }
  }
}
  
void MaxSolution::remapLabels(int max_cluster_size = 10)
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  char buf[1024];
  for (TripleVectorNonConstIt it = _triples.begin(); it != _triples.end(); ++it)
  {
    Triple& sol = *it;
    
    int n = sol._F.getNrCols();
    for (int j = 0; j < n; ++j)
    {
      StringVector s;
      boost::split(s, sol._F.getColLabel(j), boost::is_any_of(";"));
      
      std::string new_label = "";
      bool first = true;
      int count = 0;
      for (StringVectorIt it2 = s.begin(); it2 != s.end()
           && (count < max_cluster_size || max_cluster_size == -1); ++it2, ++count)
      {
        if (*it2 != "")
        {
          int org_j = boost::lexical_cast<int>(*it2);
          if (first)
          {
            first = false;
          }
          else
          {
            new_label += "\\n";
          }
          StringVector s2 ;
          boost::split(s2, _F.getColLabel(org_j), boost::is_any_of(","));
          new_label += s2.front();
        }
      }
      if (count == max_cluster_size && s.size() > max_cluster_size)
      {
        snprintf(buf, 1024, "\\n[%d more]", (int)s.size() - max_cluster_size);
        new_label += buf;
      }
      sol._F.setColLabel(j, new_label);
    }
  }
}
  
std::ostream& operator<<(std::ostream& out,
                         const MaxSolution& solution)
{
  int size = solution.size();
  out << size << " #sols" << std::endl;
  out << std::endl;

  // input matrix
  out << solution.getF();
  out << std::endl;
  
  for (int i = 0; i < size; ++i)
  {
    out << solution.solution(i);
    out << std::endl;
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         MaxSolution& solution)
{
  std::string line;
  
  vaff::getline(in, line);
  std::stringstream ss(line);
  
  int size;
  ss >> size;
  
  vaff::getline(in, line);
  
  in >> solution._F;
  solution._F.setLabels(in);
  
  vaff::getline(in, line);
  
  for (int i = 0; i < size; ++i)
  {
    MaxSolution::Triple sol;
    in >> sol;
    solution.add(sol);
    vaff::getline(in, line);
  }
  
  return in;
}
  
} // namespace vaff
