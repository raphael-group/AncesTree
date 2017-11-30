/*
 *  solutiongraph.h
 *
 *   Created on: 4-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef SOLUTIONGRAPH_H
#define SOLUTIONGRAPH_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include "utils.h"
#include "realmatrix.h"
#include "maxsolution.h"
#include "clonaltree.h"

namespace vaff {
  
class SolutionGraph
{
public:
  typedef lemon::ListBpGraph MixingGraph;
  typedef ClonalTree::Tree Tree;
  DIGRAPH_TYPEDEFS(ClonalTree::Tree);
  typedef MixingGraph::Node BpNode;
  typedef MixingGraph::NodeIt BpNodeIt;
  typedef MixingGraph::BlueNode BpBlueNode;
  typedef MixingGraph::BlueNodeIt BpBlueNodeIt;
  typedef MixingGraph::RedNode BpRedNode;
  typedef MixingGraph::RedNodeIt BpRedNodeIt;
  typedef MixingGraph::Edge BpEdge;
  typedef MixingGraph::EdgeIt BpEdgeIt;
  typedef MixingGraph::RedNodeMap<int> IntBpRedNodeMap;
  typedef MixingGraph::BlueNodeMap<int> IntBpBlueNodeMap;
  typedef MixingGraph::NodeMap<StlBoolVector> BoolVectorBpNodeMap;
  
  SolutionGraph(const MaxSolution::Triple& sol,
                double threshold,
                double beta);
  
  void writeDOT(std::ostream& out) const;
  
  const MixingGraph& getMixingGraph() const
  {
    return _G;
  }
  
  void writeEdgeList(std::ostream& out) const;
  
  void writeLeaves(std::ostream& out) const;
    
private:
  typedef std::vector<Node> NodeVector;
  typedef Tree::NodeMap<StlBoolVector> BoolVectorNodeMap;
  typedef Tree::NodeMap<BpBlueNode> BpBlueNodeNodeMap;
  
  typedef std::vector<BpRedNode> BpRedNodeVector;
  typedef std::vector<BpBlueNode> BpBlueNodeVector;
  typedef MixingGraph::RedNodeMap<StlBoolVector> BoolVectorRedNodeMap;
  typedef MixingGraph::BlueNodeMap<Node> NodeBpBlueNodeMap;
  typedef MixingGraph::EdgeMap<double> MixingEdgeMap;
  
  const MaxSolution::Triple& _sol;

  BpBlueNodeNodeMap _toMixingGraph;
  BoolNodeMap _leaf;
  
  // red node: samples
  // blue node: deconvoluted samples
  lemon::ListBpGraph _G;
  IntBpRedNodeMap _bpNodeToRow;
  BpRedNodeVector _rowToBpNode;
  IntBpBlueNodeMap _bpNodeToBasisRow;
  BoolVectorBpNodeMap _mixingLabel;
  MixingEdgeMap _mixingFraction;
  NodeBpBlueNodeMap _toTree;
  double _threshold;
  double _beta;
  
  void constructMixingGraph();
};

  
} // namespace vaff

#endif // SOLUTIONGRAPH_H
