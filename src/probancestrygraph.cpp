/*
 *  probancestrygraph.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "probancestrygraph.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <lemon/bfs.h>
#include <iomanip>

namespace vaff {
  
ProbAncestryGraph::ProbAncestryGraph()
  : BaseAncestryGraph()
  , _intermediateArc(_G)
{
}
  
ProbAncestryGraph::ProbAncestryGraph(const AncestryMatrix& A,
                                     const ReadCountMatrix& R,
                                     double alpha,
                                     double gamma)
  : BaseAncestryGraph(A.getNrCols())
  , _intermediateArc(_G)
{
  constructGraph(A, R, alpha, gamma);
}
  
void ProbAncestryGraph::constructGraph(const AncestryMatrix& A,
                                       const ReadCountMatrix& R,
                                       double alpha,
                                       double gamma)
{
  RealIntervalMatrix CI;
  R.computeConfidenceIntervals(CI, gamma);
  
  int n = 0;
  int nn = A.getNrCols();
  for (int j = 0; j < nn; ++j)
  {
    bool feasible = true;
    for (int i = 0; i < CI.getNrRows(); ++i)
    {
      if (CI(i,j).first > 0.5)
      {
        feasible = false;
        break;
      }
    }
    
    Node v_j = lemon::INVALID;
    if (feasible)
    {
      v_j = _G.addNode();
      _nodeToColumn[v_j] = j;
      ++n;
    }
    _columnToNode[j] = v_j;
  }
  
  for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
  {
    int j = _nodeToColumn[v_j];
    for (NodeIt v_k(_G); v_k != lemon::INVALID; ++v_k)
    {
      int k = _nodeToColumn[v_k];
      if (j == k)
      {
        continue;
      }
      
      double prob_j_precedes_k = A(j, k);
      bool j_precedes_k = prob_j_precedes_k >= 0.5 - alpha;
      if (j_precedes_k)
      {
        Arc a = _G.addArc(v_j, v_k);
        _intermediateArc[a] = prob_j_precedes_k <= 0.5 + alpha;
        _prob[a] = prob_j_precedes_k;
      }
    }
  }
}
  
void ProbAncestryGraph::removeCycles(const AncestryMatrix& A,
                                     double alpha,
                                     StlIntMatrix& toOrgColumns) const
{
  // hide arcs with probability > 0.5 + alpha
  IntNodeMap sccMap(_G);
  int n = lemon::stronglyConnectedComponents(lemon::filterArcs(_G, _intermediateArc), sccMap);
  
  toOrgColumns = StlIntMatrix(n, StlIntVector());
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    toOrgColumns[sccMap[v]].push_back(_nodeToColumn[v]);
  }
}
  
void ProbAncestryGraph::contract(const AncestryMatrix& A,
                                 const StlIntMatrix& toOrginalColumns,
                                 double beta,
                                 ProbAncestryGraph& H) const
{
  H._G.clear();

  int n = toOrginalColumns.size();
  H._G.reserveNode(n);
  H._columnToNode = NodeVector(n, lemon::INVALID);
  
  IntNodeMap toNewColumn(_G);
  for (int i = 0; i < n; ++i)
  {
    Node v_i = H._G.addNode();
    H._nodeToColumn[v_i] = i;
    H._columnToNode[i] = v_i;
    
    const StlIntVector& S = toOrginalColumns[i];
    for (StlIntVectorIt it = S.begin(); it != S.end(); ++it)
    {
      toNewColumn[_columnToNode[*it]] = i;
    }
  }
  
  StlBoolMatrix bitmap(n, StlBoolVector(n, false));
  for (int i = 0; i < n; ++i)
  {
    bitmap[i][i] = true;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    
    int index_u = _nodeToColumn[u];
    int index_v = _nodeToColumn[v];
    
    int index_uu = toNewColumn[u];
    int index_vv = toNewColumn[v];
    
    Node uu = H._columnToNode[index_uu];
    Node vv = H._columnToNode[index_vv];
    
    if (index_uu != index_vv
        && !bitmap[index_uu][index_vv]
        && A(index_u, index_v) >= beta)
    {
      bitmap[index_uu][index_vv] = true;
      Arc a = H._G.addArc(uu, vv);
      H._prob[a] = A(index_u, index_v);
    }
  }
  
//  for (ArcIt a(H._G); a != lemon::INVALID; ++a)
//  {
//    std::cout << H._nodeToColumn[H._G.source(a)] << " -> " << H._nodeToColumn[H._G.target(a)] << std::endl;
//  }
}
  
int ProbAncestryGraph::numberOfNodesInfCI() const
{
  int res = 0;

  int n = _columnToNode.size();
  for (int j = 0; j < n; ++j)
  {
    if (_columnToNode[j] == lemon::INVALID)
      ++res;
  }
  
  return res;
}
  
void ProbAncestryGraph::writeDOT(const ReadCountMatrix& R,
                                 const StlIntVector& S,
                                 std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  out << "\tlayout=neato" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << std::setprecision(2);
  
  // only show intra cluster edges
  lemon::ArcLookUp<Digraph> arcLookup(_G);
  
  for (StlIntVectorIt it1 = S.begin(); it1 != S.end(); ++it1)
  {
    int j = *it1;
    out << "\t" << j << " [label=\"" << R.getRowLabel(j) << "\"]" << std::endl;
    
    for (StlIntVectorIt it2 = S.begin(); it2 != S.end(); ++it2)
    {
      if (it1 == it2) continue;
      int k = *it2;
      Arc a = arcLookup(_columnToNode[j], _columnToNode[k]);
      
      if (a != lemon::INVALID)
      {
        out << "\t" << j << " -> " << k << " " << " ["
            << (_intermediateArc[a] ? "color=red," : "")
            << "minlen=10,penwidth="
            << 1 + 1 * (fabs(0.5 - _prob[a])) << ",label=\"" << _prob[a] << "\"]" << std::endl;
      }
    }
  }
  
  out << "}" << std::endl;
}
  
void ProbAncestryGraph::writeDOT(const ReadCountMatrix& R,
                                 const StlIntMatrix& M,
                                 std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << "\tlayout=neato" << std::endl;
  
  int color = 0;
  for (StlIntMatrixIt it = M.begin(); it != M.end(); ++it)
  {
    const StlIntVector& S = *it;
    bool fill = S.size() > 1;
    for (StlIntVectorIt it2 = S.begin(); it2 != S.end(); ++it2)
    {
      int j = *it2;
      out << "\t" << j << " [label=\"" << R.getRowLabel(j) << "\"";
      if (fill)
      {
        out << ",colorscheme=paired10,style=filled,color=" << color + 1;
      }
      out << "]" << std::endl;
    }
    if (fill)
    {
      color = (color + 1) % 10;
    }
  }
  
  for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
  {
    int j = _nodeToColumn[v_j];
    out << "\t" << j << " [label=\"" << R.getRowLabel(j) << "\"";
    out << "]" << std::endl;
  }
  
  lemon::ArcLookUp<Digraph> arcLookup(_G);
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_j = _G.source(a);
    int j = _nodeToColumn[v_j];
    Node v_k = _G.target(a);
    int k = _nodeToColumn[v_k];
    
    out << "\t" << j << " -> " << k << " [label=" << std::setprecision(2) << _prob[a]
        << ",penwidth=" << 5 * _prob[a];
    Arc rev_a = arcLookup(v_k, v_j);
    if (_intermediateArc[a] && rev_a != lemon::INVALID && _intermediateArc[rev_a])
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }
  out << "}" << std::endl;
}
  
void ProbAncestryGraph::writeDOT(const ReadCountMatrix& R,
                                 std::ostream& out) const
{
//  BoolNodeMap infeasible(_G, false);
//  
  out << "digraph G {" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << "\tlayout=neato" << std::endl;
  for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
  {
    int j = _nodeToColumn[v_j];
    out << "\t" << j << " [label=\"" << R.getRowLabel(j) << "\"";
    out << "]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_j = _G.source(a);
    int j = _nodeToColumn[v_j];
    Node v_k = _G.target(a);
    int k = _nodeToColumn[v_k];
    
    out << "\t" << j << " -> " << k << " [label=" << _prob[a];
    if (_intermediateArc[a])
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }
  out << "}" << std::endl;
}
  
void ProbAncestryGraph::createLabels(const StlIntMatrix& toOrginalColumns,
                                     StringVector& label)
{
  const int n = toOrginalColumns.size();
  char buf[1024];
  
  label = StringVector(n, "");
  
  for (int i = 0; i < n; ++i)
  {
    const StlIntVector& S = toOrginalColumns[i];
    std::string lbl;
    bool first = true;
    for (StlIntVectorIt it = S.begin(); it != S.end(); ++it)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        lbl += ";";
      }

      snprintf(buf, 1024, "%d", *it);
      lbl += buf;
    }
    label[i] = lbl;
  }
}
  
void ProbAncestryGraph::printClusterSize(const StlIntMatrix& toOrginalColumns,
                                         std::ostream& out)
{
  const int n = toOrginalColumns.size();
  for (int i = 0; i < n; ++i)
  {
    out << toOrginalColumns[i].size() << std::endl;
  }
}
  
void ProbAncestryGraph::printInterClusterCoherence(const Digraph& G,
                                                   const IntNodeMap& nodeToMutation,
                                                   const StringVector& mutationToLabel,
                                                   const AncestryMatrix& A,
                                                   std::ostream& out)
{
  for (ArcIt a(G); a != lemon::INVALID; ++a)
  {
    Node v_j = G.source(a);
    Node v_k = G.target(a);
    
    int j = nodeToMutation[v_j];
    int k = nodeToMutation[v_k];
    
    const std::string& label_v_j = mutationToLabel[j];
    const std::string& label_v_k = mutationToLabel[k];
    
    StringVector s_j, s_k;
    boost::split(s_j, label_v_j, boost::is_any_of(";"));
    boost::split(s_k, label_v_k, boost::is_any_of(";"));
    
    double p = 0;
    for (StringVectorIt it1 = s_j.begin(); it1 != s_j.end(); ++it1)
    {
      int org_j = boost::lexical_cast<int>(*it1);
      for (StringVectorIt it2 = s_k.begin(); it2 != s_k.end(); ++it2)
      {
        int org_k = boost::lexical_cast<int>(*it2);
        p = std::max(p, A(org_j, org_k));
      }
    }
    
    out << p << std::endl;
  }
}
  
double ProbAncestryGraph::interClusterCoherence(const Digraph& G,
                                                const IntNodeMap& nodeToMutation,
                                                const StringVector& mutationToLabel,
                                                const AncestryMatrix& A)
{
  int total_count = 0;
  int anc_count = 0;
  for (ArcIt a(G); a != lemon::INVALID; ++a)
  {
    Node v_j = G.source(a);
    Node v_k = G.target(a);
    
    int j = nodeToMutation[v_j];
    int k = nodeToMutation[v_k];
    
    const std::string& label_v_j = mutationToLabel[j];
    const std::string& label_v_k = mutationToLabel[k];
    
    StringVector s_j, s_k;
    boost::split(s_j, label_v_j, boost::is_any_of(";"));
    boost::split(s_k, label_v_k, boost::is_any_of(";"));
    
    double p = 0;
    for (StringVectorIt it1 = s_j.begin(); it1 != s_j.end(); ++it1)
    {
      int org_j = boost::lexical_cast<int>(*it1);
      for (StringVectorIt it2 = s_k.begin(); it2 != s_k.end(); ++it2)
      {
        ++total_count;
        int org_k = boost::lexical_cast<int>(*it2);
        p = std::max(p, A(org_j, org_k));
      }
    }
    
    if (p >= 0.7)
    {
      return ++anc_count;
    }
  }
  
  return (double) anc_count / (double) total_count;
}
  
void ProbAncestryGraph::printIntraClusterCoherence(const Digraph& G,
                                                   const IntNodeMap& nodeToMutation,
                                                   const StringVector& mutationToLabel,
                                                   const AncestryMatrix& A,
                                                   std::ostream& out)
{
  for (NodeIt v_j(G); v_j != lemon::INVALID; ++v_j)
  {
    int j = nodeToMutation[v_j];
    std::string label_v_j = mutationToLabel[j];
    
    StringVector s;
    boost::split(s, label_v_j, boost::is_any_of(";"));
    
    for (StringVectorIt it1 = s.begin(); it1 != s.end(); ++it1)
    {
      int org_j_1 = boost::lexical_cast<int>(*it1);
      for (StringVectorIt it2 = it1 + 1; it2 != s.end(); ++it2)
      {
        int org_j_2 = boost::lexical_cast<int>(*it2);
        double p = std::max(A(org_j_1, org_j_2), A(org_j_2, org_j_1));
        out << p << " " << label_v_j << " " << *it1 << " " << *it2 << std::endl;
//        out << A(org_j_1, org_j_2) << " " << label_v_j << " " << *it1 << " " << *it2 << std::endl;
//        out << A(org_j_2, org_j_1) << " " << label_v_j << " " << *it2 << " " << *it1 << std::endl;
      }
    }
  }
}
  
void ProbAncestryGraph::intraClusterCoherence(const Digraph& G,
                                              const IntNodeMap& nodeToMutation,
                                              const StringVector& mutationToLabel,
                                              const AncestryMatrix& A,
                                              double& incomparable,
                                              double& ancestral,
                                              double& comparable)
{
  int total_count = 0;
  int inc_count = 0;
  int anc_count = 0;
  int com_count = 0;
  for (NodeIt v_j(G); v_j != lemon::INVALID; ++v_j)
  {
    int j = nodeToMutation[v_j];
    std::string label_v_j = mutationToLabel[j];
    
    StringVector s;
    boost::split(s, label_v_j, boost::is_any_of(";"));
    
    for (StringVectorIt it1 = s.begin(); it1 != s.end(); ++it1)
    {
      int org_j_1 = boost::lexical_cast<int>(*it1);
      for (StringVectorIt it2 = it1 + 1; it2 != s.end(); ++it2)
      {
        ++total_count;
        int org_j_2 = boost::lexical_cast<int>(*it2);
        double p = std::max(A(org_j_1, org_j_2), A(org_j_2, org_j_1));
        if (p <= 0.05)
        {
          ++inc_count;
        }
        else if (p >= 0.7)
        {
          ++anc_count;
        }
        else
        {
          ++com_count;
        }
      }
    }
  }
  
  if (total_count == 0)
  {
    incomparable = ancestral = 0;
    comparable = 1;
  }
  else
  {
    incomparable = (double) inc_count / (double) total_count;
    ancestral = (double) anc_count / (double) total_count;
    comparable = (double) com_count / (double) total_count;
  }
}
  
} // namespace vaff
