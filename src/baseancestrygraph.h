/*
 *  baseancestrygraph.h
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef BASEANCESTRYGRAPH_H
#define BASEANCESTRYGRAPH_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <lemon/connectivity.h>
#include "utils.h"
#include "realmatrix.h"
#include "ppmatrix.h"

namespace vaff {
  
class BaseAncestryGraph
{
public:
  typedef lemon::ListDigraph Digraph;
  DIGRAPH_TYPEDEFS(Digraph);
  typedef std::vector<Node> NodeVector;
  
  BaseAncestryGraph();
  
  BaseAncestryGraph(int n);
  
  void writeDOT(std::ostream& out) const;
  
  bool isTransitive() const
  {
    lemon::ArcLookUp<Digraph> arcLookUp(_G);
    for (ArcIt a(_G); a != lemon::INVALID; ++a)
    {
      Node v_i = _G.source(a);
      Node v_j = _G.target(a);
      for (OutArcIt aa(_G, v_j); aa != lemon::INVALID; ++aa)
      {
        Node v_k = _G.target(aa);
        if (v_k == v_i)
        {
          continue;
        }
        if (arcLookUp(v_i, v_k) == lemon::INVALID)
        {
          std::cout << _nodeToColumn[v_i] << " -> " << _nodeToColumn[v_j] << " -> " << _nodeToColumn[v_k] << std::endl;
          return false;
        }
      }
    }
    return true;
  }
  
  bool isDAG() const
  {
//    IntNodeMap comp(_G);
//    int n = lemon::stronglyConnectedComponents(_G, comp);
//    
//    for (int j = 0; j < n; ++j)
//    {
//      std::cout << j << " :";
//      for (NodeIt v(_G); v != lemon::INVALID; ++v)
//      {
//        if (comp[v] == j)
//        {
//          std::cout << " " << _nodeToColumn[v];
//        }
//      }
//      std::cout << std::endl;
//    }
//    
//    std::cout << n << " / " << lemon::countNodes(_G) << std::endl;
//    std::cout << lemon::dag(_G) << std::endl;
//    std::cout << lemon::loopFree(_G) << std::endl;
    return lemon::dag(_G);
  }
  
  int numberOfNonTrivialSCC() const
  {
    IntNodeMap comp(_G);
    int nComp = lemon::stronglyConnectedComponents(_G, comp);
    
    int res = 0;
    
    for (int j = 0; j < nComp; ++j)
    {
      int res2 = 0;
      for (NodeIt v(_G); v != lemon::INVALID; ++v)
      {
        if (comp[v] == j)
        {
          ++res2;
        }
        if (res2 == 2)
        {
          ++res;
          break;
        }
      }
    }
    
    return res;
  }
  
  int numberOfIsolatedNodes() const
  {
    int res = 0;
    for (NodeIt v(_G); v != lemon::INVALID; ++v)
    {
      if (InArcIt(_G, v) == lemon::INVALID && OutArcIt(_G, v) == lemon::INVALID)
      {
        ++res;
      }
    }
    return res;
  }
  
  int numberOfNodesInDeg0() const
  {
    int res = 0;
    for (NodeIt v(_G); v != lemon::INVALID; ++v)
    {
      if (InArcIt(_G, v) == lemon::INVALID)
      {
        ++res;
      }
    }
    return res;
  }
  
  int largestArborescence() const
  {
    lemon::Bfs<Digraph> bfs(_G);
    
    int res = 0;
    for (NodeIt v(_G); v != lemon::INVALID; ++v)
    {
      if (InArcIt(_G, v) == lemon::INVALID)
      {
        bfs.init();
        bfs.run(v);
        int s = 0;
        for (NodeIt u(_G); u != lemon::INVALID; ++u)
        {
          if (bfs.reached(u))
            ++s;
        }
        res = std::max(res, s);
      }
    }
    return res;
  }
  
  int maxOutDegree() const
  {
    int res = 0;
    for (NodeIt v(_G); v != lemon::INVALID; ++v)
    {
      int d = 0;
      for (OutArcIt a(_G, v); a != lemon::INVALID; ++a, ++d);
      res = std::max(res, d);
    }
    return res;
  }
  
  const Digraph& getG() const
  {
    return _G;
  }
  
  Node mapColumnToNode(int p) const
  {
    assert(0 <= p && p < _columnToNode.size());
    return _columnToNode[p];
  }
  
  int mapNodeToColumn(Node v) const
  {
    return _nodeToColumn[v];
  }
  
  const IntNodeMap& getNodeToColumnMap() const
  {
    return _nodeToColumn;
  }
  
  const NodeVector& getColumnToNodeVector() const
  {
    return _columnToNode;
  }
  
  void contract(const StlIntMatrix& toOrginalColumns,
                BaseAncestryGraph& H) const;
  
  const DoubleArcMap& getProbMap() const
  {
    return _prob;
  }
  
  double getProb(Arc a) const
  {
    return _prob[a];
  }
  
protected:
  typedef lemon::OutDegMap<Digraph> OutDegMap;
  Digraph _G;
  NodeVector _columnToNode;
  IntNodeMap _nodeToColumn;
  OutDegMap _outDegree;
  DoubleArcMap _prob;
};
  
} // namespace vaff

#endif // BASEANCESTRYGRAPH_H
