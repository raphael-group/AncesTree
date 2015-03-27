/*
 *  clonaltree.cpp
 *
 *   Created on: 5-jan-2015
 *       Author: M. El-Kebir
 */

#include "clonaltree.h"
#include <lemon/connectivity.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace vaff {
  
ClonalTree::ClonalTree()
  : _T()
  , _root(lemon::INVALID)
  , _nodeToMutation(_T)
  , _mutationToNode()
  , _prob(_T)
  , _B()
{  
}
  
ClonalTree::ClonalTree(const ClonalTree& CT)
  : _T()
  , _root(lemon::INVALID)
  , _nodeToMutation(_T)
  , _mutationToNode(lemon::countNodes(CT._T), lemon::INVALID)
  , _prob(_T)
  , _B(CT._B)
{
  Tree::NodeMap<Node> m(CT._T);
  lemon::digraphCopy(CT._T, _T).nodeMap(CT._nodeToMutation, _nodeToMutation).nodeRef(m).arcMap(CT._prob, _prob).run();

  _root = m[CT._root];
  for (NodeIt v(CT._T); v != lemon::INVALID; ++v)
  {
    Node new_v = m[v];
    int mut_v = CT._nodeToMutation[v];
    _mutationToNode[mut_v] = new_v;
    _nodeToMutation[new_v] = mut_v;
  }
}
 
ClonalTree::ClonalTree(const PpMatrix& B)
  : _T()
  , _root(lemon::INVALID)
  , _nodeToMutation(_T)
  , _mutationToNode()
  , _prob(_T)
  , _B(B)
{
  assert(B.getNrCols() == B.getNrRows());
  constructT();
}
  
void ClonalTree::constructB(Node v)
{
  int n = _mutationToNode.size();
  int mutation_v = _nodeToMutation[v];
  
  if (v != _root)
  {
    Node u = _T.source(InArcIt(_T, v));
    int mutation_u = _nodeToMutation[u];
    for (int j = 0; j < n; ++j)
    {
      _B.set(mutation_v, j, _B(mutation_u, j));
    }
  }
  _B.set(mutation_v, mutation_v, true);
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node v = _T.target(a);
    constructB(v);
    _B.setProb(_nodeToMutation[v], _prob[a]);
  }
}

void ClonalTree::constructT()
{
  int n = _B.getNrCols();
  _root = lemon::INVALID;
  _mutationToNode = NodeVector(n, lemon::INVALID);
  
  // construct nodes and identify root node
  _T.reserveNode(n);
  _T.reserveArc(n - 1);
  for (int j = 0; j < n; ++j)
  {
    Node v = _T.addNode();
    _nodeToMutation[v] = j;
    _mutationToNode[j] = v;
    
    if (_B.rowSum(j) == 1)
    {
      assert(_root == lemon::INVALID);
      _root = _mutationToNode[j];
    }
  }
  
  // construct arcs
  for (int j1 = 0; j1 < n; ++j1)
  {
    for (int j2 = 0; j2 < n; ++j2)
    {
      if (j1 == j2)
        continue;
      
      if (_B.isRowContained(j1, j2) == 1)
      {
        Arc a = _T.addArc(_mutationToNode[j1], _mutationToNode[j2]);
        _prob[a] = _B.prob(j2);
      }
    }
  }
//  
//  assert(lemon::dag(_T));
//  std::cout << lemon::countArcs(_T) << std::endl;
}
  
RealMatrix ClonalTree::getU(const RealMatrix& F) const
{
  assert(F.getNrCols() == _B.getNrCols());

  int m = F.getNrRows();
  int n = F.getNrCols();
         
  RealMatrix U(m, n);
  
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      double u_ij = 2 * F(i, j);
      
      Node v_j = _mutationToNode[j];
      for (OutArcIt a(_T, v_j); a != lemon::INVALID; ++a)
      {
        Node v_k = _T.target(a);
        int k = _nodeToMutation[v_k];
        
        u_ij -= 2 * F(i, k);
      }
      
      if (!g_tol.nonZero(u_ij))
        u_ij = 0;
      
      U.set(i, j, u_ij);
    }
  }
  
  return U;
}

std::ostream& operator<<(std::ostream& out,
                         const ClonalTree& T)
{
  out << T._B;
  return out;
}

std::istream& operator>>(std::istream& in,
                         ClonalTree& T)
{
  in >> T._B;
  T.constructT();
  return in;
}
  
} // namespace vaff
