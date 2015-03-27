/*
 *  clonaltree.h
 *
 *   Created on: 5-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef CLONALTREE_H
#define CLONALTREE_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <vector>
#include "ppmatrix.h"
#include "realmatrix.h"
#include "utils.h"
#include "ancestrymatrix.h"

namespace vaff {
  
class ClonalTree
{
public:
  typedef lemon::ListDigraph Tree;
  
  DIGRAPH_TYPEDEFS(Tree);
  
  ClonalTree();
  

  ClonalTree(const ClonalTree& CT);
  
  template<typename TREE, typename INM, typename DAM>
  ClonalTree(const TREE& T,
             const INM& nodeToMutation,
             const DAM& prob);
  
  ClonalTree(const PpMatrix& B);
  
  ClonalTree& operator=(const ClonalTree& other)
  {
    if (this == &other)
    {
      return *this;
    }
    
    _B = other._B;
    _mutationToNode = NodeVector(_B.getNrCols(), lemon::INVALID);
    
    Tree::NodeMap<Node> m(other._T);
    lemon::digraphCopy(other._T, _T).arcMap(other._prob, _prob).nodeMap(other._nodeToMutation, _nodeToMutation).nodeRef(m).run();
    
    _root = m[other._root];
    for (NodeIt v(other._T); v != lemon::INVALID; ++v)
    {
      Node new_v = m[v];
      int mut_v = other._nodeToMutation[v];
      _mutationToNode[mut_v] = new_v;
      _nodeToMutation[new_v] = mut_v;
    }
    
    return *this;
  }
  
  bool operator==(const ClonalTree& other) const
  {
    return _B == other._B;
  }
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const ClonalTree& T);
  
  friend std::istream& operator>>(std::istream& in,
                                  ClonalTree& T);
  
  const Tree& getT() const
  {
    return _T;
  }
  
  const PpMatrix& getB() const
  {
    return _B;
  }

  bool ancestral(int j, int k) const
  {
    return j != k && _B(k, j);
  }
  
  bool ancestral(Node v_j, Node v_k) const
  {
    int j = _nodeToMutation[v_j];
    int k = _nodeToMutation[v_k];
    
    return ancestral(j, k);
  }
  
  bool incomparable(int j, int k) const
  {
    return !_B(j, k) && !_B(k, j);
  }
  
  bool incomparable(Node v_j, Node v_k) const
  {
    int j = _nodeToMutation[v_j];
    int k = _nodeToMutation[v_k];
    
    return incomparable(j, k);
  }
  
  Node getRoot() const
  {
    return _root;
  }
  
  Node getParent(Node v) const
  {
    InArcIt a(_T, v);
    if (a == lemon::INVALID)
    {
      return lemon::INVALID;
    }
    return _T.source(a);
  }
  
  Node mutationToNode(int j) const
  {
    assert(0 <= j && j < _mutationToNode.size());
    return _mutationToNode[j];
  }
  
  int nodeToMutation(Node v) const
  {
    return _nodeToMutation[v];
  }
  
  const IntNodeMap& getNodeToMutationMap() const
  {
    return _nodeToMutation;
  }
  
  const DoubleArcMap& getProbMap() const
  {
    return _prob;
  }
  
  RealMatrix getU(const RealMatrix& F) const;
  
private:
  typedef std::vector<Node> NodeVector;
  
  Tree _T;
  Node _root;
  IntNodeMap _nodeToMutation;
  NodeVector _mutationToNode;
  DoubleArcMap _prob;
  
  PpMatrix _B;
  
  void constructB(Node v);
  
  void constructT();
};

template<typename TREE, typename INM, typename DAM>
ClonalTree::ClonalTree(const TREE& T,
                       const INM& nodeToMutation,
                       const DAM& prob)
  : _T()
  , _root(lemon::INVALID)
  , _nodeToMutation(_T)
  , _mutationToNode(lemon::countNodes(T), lemon::INVALID)
  , _prob(_T)
  , _B(_mutationToNode.size(), _mutationToNode.size())
{
  lemon::digraphCopy(T, _T).arcMap(prob, _prob).nodeMap(nodeToMutation, _nodeToMutation).run();
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    if (InArcIt(_T, v) == lemon::INVALID)
    {
      _root = v;
    }
    _mutationToNode[_nodeToMutation[v]] = v;
  }

  constructB(_root);
}
  
} // namespace vaff

#endif // CLONALTREE_H
