/*
 *  solutiongraph.cpp
 *
 *   Created on: 4-jan-2015
 *       Author: M. El-Kebir
 */

#include "solutiongraph.h"
#include <lemon/bfs.h>
#include <iomanip>

namespace vaff {

SolutionGraph::SolutionGraph(const MaxSolution::Triple& sol,
                             double threshold,
                             double beta)
  : _sol(sol)
  , _toMixingGraph(sol._T.getT())
  , _leaf(sol._T.getT())
  , _G()
  , _bpNodeToRow(_G)
  , _rowToBpNode()
  , _bpNodeToBasisRow(_G)
  , _mixingLabel(_G)
  , _mixingFraction(_G)
  , _toTree(_G)
  , _threshold(threshold)
  , _beta(beta)
{
  for (NodeIt v(sol._T.getT()); v != lemon::INVALID; ++v)
  {
    _leaf[v] = OutArcIt(sol._T.getT(), v) == lemon::INVALID;
  }
  constructMixingGraph();
}

void SolutionGraph::constructMixingGraph()
{
  const Tree& T = _sol._T.getT();
  int m = _sol._F.getNrRows();
  int n = _sol._F.getNrCols();
  lemon::mapFill(T, _toMixingGraph, lemon::INVALID);
  
  // let's first add the sample nodes
  _rowToBpNode = BpRedNodeVector(m, lemon::INVALID);
  for (int i = 0; i < m; ++i)
  {
    BpRedNode v = _G.addRedNode();
    _rowToBpNode[i] = v;
    _bpNodeToRow[v] = i;
    _mixingLabel[v] = discretize(_sol._F.getMatrix()[i]);
  }
  
  // now let's add the mixing nodes
  int mm = 0;
  for (int j = 0; j < n; ++j)
  {
    for (int i = 0; i < m; ++i)
    {
      if (_sol._U(i, j) != 0 && _sol._U(i, j) >= _threshold)
      {
        ++mm;
        
        BpBlueNode v = _G.addBlueNode();
        _bpNodeToBasisRow[v] = j;
        _mixingLabel[v] = _sol._T.getB().getMatrix()[j];
        bool found = false;
        for (NodeIt vv(T); vv != lemon::INVALID; ++vv)
        {
          if (_sol._T.getB().getMatrix()[_sol._T.nodeToMutation(vv)] == _mixingLabel[v])
          {
            assert(!found);
            _toTree[v] = vv;
            _toMixingGraph[vv] = v;
//            std::cerr << _sol._F.getColLabel(j) << std::endl;
            found = true;
          }
        }
        assert(found);
        
        break;
      }
    }
  }
  
//  for (int i = 0; i < mm; ++i)
//  {
//    BpBlueNode v = _G.addBlueNode();
//    _bpNodeToBasisRow[v] = i;
//    _mixingLabel[v] = _sol._T.getB().getMatrix()[i];
//    bool found = false;
//    for (NodeIt vv(T); vv != lemon::INVALID; ++vv)
//    {
//      if (_sol._T.getB().getMatrix()[_sol._T.nodeToMutation(vv)] == _mixingLabel[v])
//      {
//        assert(!found);
//        _toTree[v] = vv;
//        _toMixingGraph[vv] = v;
//        std::cerr << _sol._F.getColLabel(_sol._T.nodeToMutation(vv)) << std::endl;
//        found = true;
//      }
//    }
//    assert(found);
//  }
  
  // and finally the mixing edges
  int mmm = _sol._U.getNrRows();
  int nnn = _sol._U.getNrCols();
  for (int i = 0; i < mmm; ++i)
  {
    for (int j = 0; j < nnn; ++j)
    {
      double u_ij = _sol._U(i, j);
      if (u_ij && u_ij >= _threshold)
      {
        BpEdge e = _G.addEdge(_rowToBpNode[i], _toMixingGraph[_sol._T.mutationToNode(j)]);
        _mixingFraction[e] = u_ij;
      }
    }
  }
}

void SolutionGraph::writeDOT(std::ostream& out) const
{
  static int fontsizeBox = 60;
  static int fontsize = 35;
  static int minPenwidth = 5;
  
  const Tree& T = _sol._T.getT();
  
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "graph G {" << std::endl;
  
  out << "\tsubgraph mixed {" << std::endl;
  for (BpRedNodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t\ts" << _G.id(v)
        << " [colorscheme=paired10,penwidth=" << minPenwidth
        << ",fontsize=" << fontsizeBox
        << ",color=" << (_G.id(v) % 10) + 1 << ",shape=box,label=\"";
    out << _sol._F.getRowLabel(_bpNodeToRow[v]) << "\"]" << std::endl;
  }
  out << "\t}" << std::endl;
  
  lemon::Bfs<Tree> bfs(T);
  bfs.run(_sol._T.getRoot());
  int maxLevel = lemon::mapMaxValue(T, bfs.distMap());
  
  out << "\tsubgraph unmixed {" << std::endl;
  for (int l = 0; l < maxLevel; ++l)
  {
    out << "\t\tsubgraph " << l << " {" << std::endl;
    out << "\t\t\trank=same" << std::endl;
    for (NodeIt v(T); v != lemon::INVALID; ++v)
    {
      if (bfs.dist(v) == l && !_leaf[v])
      {
        out << "\t\t\t" << T.id(v)
            << " [penwidth=" << minPenwidth
            << ",fontsize=" << fontsize << ",label=\""
            << _sol._F.getColLabel(_sol._T.nodeToMutation(v))
            << "\"]" << std::endl;
      }
    }
    out << "\t\t}" << std::endl;
  }
  
  // leaves
  out << "\t\tsubgraph " << "leaves" << " {" << std::endl;
  out << "\t\t\trank=same" << std::endl;
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (_leaf[v])
    {
      out << "\t\t\t" << T.id(v)
          << " [penwidth=" << minPenwidth
          << ",fontsize=" << fontsize
          << ",label=\"" << _sol._F.getColLabel(_sol._T.nodeToMutation(v))
          << "\"]" << std::endl;
    }
  }
  for (BpBlueNodeIt v(_G); v != lemon::INVALID; ++v)
  {
    Node vv = _toTree[v];
    if (vv != lemon::INVALID && !_leaf[vv])
    {
      assert(_toMixingGraph[vv] == v);
      out << "\t\t\tdup" << _G.id(v)
          << " [penwidth=" << minPenwidth
          << ",fontsize=" << fontsize
          << ",label=\"" << _sol._F.getColLabel(_sol._T.nodeToMutation(vv))
          << "\"]" << std::endl;
    }
  }
  out << "\t\t}" << std::endl;
  out << "\t}" << std::endl;
  
  for (ArcIt a(T); a != lemon::INVALID; ++a)
  {
    double p_a = _sol._T.getProbMap()[a];
    
    Node u = T.source(a);
    Node v = T.target(a);
    out << "\t" << T.id(u) << " -- " << T.id(v)
        << " [penwidth=" << minPenwidth + minPenwidth * (p_a - _beta) / (1 - _beta)
        << ",fontsize=" << fontsize
        << ",label=\" " << std::setprecision(2) << p_a << "\"]" << std::endl;
  }
  
  for (BpEdgeIt e(_G); e != lemon::INVALID; ++e)
  {
    assert(_G.valid(e));
    if (_mixingFraction[e] < _threshold)
      continue;
    
    BpBlueNode v = _G.blueNode(e);
    BpRedNode u = _G.redNode(e);
    Node vv = _toTree[v];
    if (vv != lemon::INVALID && _leaf[_toTree[v]])
    {
      out << "\t" << T.id(vv) << " -- "
          << "s" << _G.id(u);
    }
    else
    {
      out << "\tdup" << _G.id(v) << " -- "
          << "s" << _G.id(u);
    }
    out << " [splines=none,colorscheme=paired10,color=" << (_G.id(u) % 10) + 1 << ",minlen=4,fontsize="
        << fontsize << ",label=" << std::setprecision(2)
        << _mixingFraction[e] << ",penwidth="
        << minPenwidth + 25 * _mixingFraction[e] << "]" << std::endl;
  }
  
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    BpBlueNode vv = _toMixingGraph[v];
    if (!_leaf[v] && vv != lemon::INVALID)
    {
      out << "\t" << T.id(v) << " -- " << "dup" << _G.id(vv)
          << " [penwidth=" << minPenwidth << ",style=dashed]" << std::endl;
    }
  }
  
  out << "}" << std::endl;
}

} // namespace vaff
