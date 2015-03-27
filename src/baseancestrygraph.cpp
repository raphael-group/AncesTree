/*
 *  baseancestrygraph.cpp
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#include "baseancestrygraph.h"
#include <lemon/bfs.h>

namespace vaff {
  
BaseAncestryGraph::BaseAncestryGraph()
  : _G()
  , _columnToNode()
  , _nodeToColumn(_G)
  , _outDegree(_G)
  , _prob(_G)
{
}
  
BaseAncestryGraph::BaseAncestryGraph(int n)
  : _G()
  , _columnToNode(n, lemon::INVALID)
  , _nodeToColumn(_G)
  , _outDegree(_G)
  , _prob(_G)
{
}
  
void BaseAncestryGraph::writeDOT(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t" << _nodeToColumn[v] << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node s = _G.source(a);
    Node t = _G.target(a);
    
    out << "\t" << _nodeToColumn[s] << " -> " << _nodeToColumn[t] << std::endl;
  }
  out << "}" << std::endl;
}
  
void BaseAncestryGraph::contract(const StlIntMatrix& toOrginalColumns,
                                 BaseAncestryGraph& H) const
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
    
    int index_uu = toNewColumn[u];
    int index_vv = toNewColumn[v];
    
    Node uu = H._columnToNode[index_uu];
    Node vv = H._columnToNode[index_vv];
    
    if (!bitmap[index_uu][index_vv])
    {
      bitmap[index_uu][index_vv] = true;
      H._G.addArc(uu, vv);
    }
  }
  
//  for (ArcIt a(H._G); a != lemon::INVALID; ++a)
//  {
//    std::cout << H._nodeToColumn[H._G.source(a)] << " -> " << H._nodeToColumn[H._G.target(a)] << std::endl;
//  }
}
  
} // namespace vaff