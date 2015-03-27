/*
 *  intmaxilpsolver.cpp
 *
 *   Created on: 8-jan-2015
 *       Author: M. El-Kebir
 */

#include "intmaxilpsolver.h"
#include <lemon/bfs.h>
#include <lemon/adaptors.h>

namespace vaff {
  
IntMaxIlpSolver::IntMaxIlpSolver(const BaseAncestryGraph& G,
                                 const RealIntervalMatrix& F_interval_clustered,
                                 const RealMatrix& F_point_unclustered,
                                 const StlIntMatrix& toUnclusteredColumn,
                                 int timeLimit)
  : _G(G)
  , _F_interval_clustered(F_interval_clustered)
  , _F_point_unclustered(F_point_unclustered)
  , _toUnclusteredColumn(toUnclusteredColumn)
  , _timeLimit(timeLimit)
  , _nodeCount(lemon::countNodes(_G.getG()))
  , _nodeToIndex(_G.getNodeToColumnMap())
  , _indexToNode(_G.getColumnToNodeVector())
  , _nodeToRootArcIndex(_G.getG(), -1)
  , _arcCount(lemon::countArcs(_G.getG()))
  , _arcToIndex(_G.getG(), -1)
  , _indexToArc(_arcCount + _nodeCount, lemon::INVALID)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _x()
  , _f()
  , _fx()
  , _g()
{
  assert(_G.isDAG());
  assert(F_interval_clustered.getNrRows() == F_point_unclustered.getNrRows());
  assert(F_interval_clustered.getNrCols() <= F_point_unclustered.getNrCols());
  
  initVariables();
  initConstraints();
  initObjective();
}
  
IntMaxIlpSolver::~IntMaxIlpSolver()
{
  _env.end();
}
  
bool IntMaxIlpSolver::solve(MaxSolution& solution)
{
  _cplex.setOut(std::cerr);
  _cplex.setWarning(std::cerr);
  _cplex.setError(std::cerr);
  
  _cplex.setParam(IloCplex::SolnPoolAGap, 0.0);
  _cplex.setParam(IloCplex::SolnPoolIntensity, 4);
  _cplex.setParam(IloCplex::PopulateLim, 20000);
  
  if (_timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, _timeLimit);
  }
  
  if (!_cplex.solve())
  {
    return false;
  }
  
  double obj_value = _cplex.getObjValue();
  
  int nSol = _cplex.getSolnPoolNsolns();
  solution.clear();
  
  for (int solIdx = 0; solIdx < nSol; ++solIdx)
  {
    if (_cplex.getObjValue(solIdx) == obj_value)
    {
      MaxSolution::Triple sol;
//      printVariables(solIdx, std::cout);
      processSolution(solIdx, sol);
      
      if (!solution.present(sol))
      {
        solution.add(sol);
      }
    }
  }
  
  std::cerr << "[" << _cplex.getObjValue() << ", " << _cplex.getBestObjValue() << "]" << std::endl;
  
  return true;
}
  
void IntMaxIlpSolver::printVariables(int solIdx, std::ostream& out) const
{
  out << "Solution " << solIdx << std::endl;
  out << "Obj value: " << _cplex.getObjValue(solIdx) << std::endl;
  
  for (int i = 0; i < _arcCount + _nodeCount; ++i)
  {
    bool val = !(fabs(_cplex.getValue(_x[i], solIdx)) <= 1e-3);
    if (val)
      out << _x[i].getName() << " = " << val << std::endl;
  }
  
  const int m = _F_interval_clustered.getNrRows();
  const int n = _F_interval_clustered.getNrCols();
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      out << _f[i][j].getName() << " = " << _cplex.getValue(_f[i][j], solIdx) << std::endl;
    }
  }
}
  
void IntMaxIlpSolver::processSolution(int solIdx, MaxSolution::Triple& sol) const
{ 
  const Digraph& G = _G.getG();
  
  // get arborescence from the solution
  StlIntVector mutationSolVec;
  IntNodeMap mutationSolMap(G, -1);
  BoolNodeMap mutationSolMap2(G, false);
  
  // find root
  for (NodeIt v(G); v != lemon::INVALID; ++v)
  {
    bool val = !(fabs(_cplex.getValue(_x[_nodeToRootArcIndex[v]], solIdx)) <= 1e-3);
    if (val)
    {
      mutationSolMap[v] = mutationSolVec.size();
      mutationSolVec.push_back(_nodeToIndex[v]);
      mutationSolMap2[v] = true;
    }
  }
  
  BoolArcMap arcSolMap(G, false);
  for (ArcIt a(G); a != lemon::INVALID; ++a)
  {
    bool val = !(fabs(_cplex.getValue(_x[_arcToIndex[a]], solIdx)) <= 1e-3);
    arcSolMap[a] = val;
    
    if (val)
    {
      Node u = G.source(a);
      if (mutationSolMap[u] == -1)
      {
        mutationSolMap[u] = mutationSolVec.size();
        mutationSolVec.push_back(_nodeToIndex[u]);
        mutationSolMap2[u] = true;
      }
      
      Node v = G.target(a);
      if (mutationSolMap[v] == -1)
      {
        mutationSolMap[v] = mutationSolVec.size();
        mutationSolVec.push_back(_nodeToIndex[v]);
        mutationSolMap2[v] = true;
      }
    }
  }
  
  const int m = _F_interval_clustered.getNrRows();
  const int n = mutationSolVec.size();
  
  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  SubDigraph subG(_G.getG(), mutationSolMap2, arcSolMap);
  
  for (SubDigraph::NodeIt v_j(subG); v_j != lemon::INVALID; ++v_j)
  {
    int j = _nodeToIndex[v_j];
    for (int i = 0; i < m; ++i)
    {
      const RealInterval& interval_ij = _F_interval_clustered(i,j);
      double f_ij = _cplex.getValue(_f[i][j]);
      assert(!g_tol.less(f_ij, interval_ij.first) && !g_tol.less(interval_ij.second, f_ij));
      double sum = 0;
      for (SubDigraph::OutArcIt a_jk(subG, v_j); a_jk != lemon::INVALID; ++a_jk)
      {
        Node v_k = G.target(a_jk);
        int k = _nodeToIndex[v_k];
        double f_ik = _cplex.getValue(_f[i][k]);
        const RealInterval& interval = _F_interval_clustered(i,k);
        assert(!g_tol.less(f_ik, interval.first) && !g_tol.less(interval.second, f_ik));
        sum += f_ik;
      }
      assert(!g_tol.less(f_ij, sum));
    }
  }
  
  sol._T = ClonalTree(subG, mutationSolMap, _G.getProbMap());

  sol._F = RealMatrix(m, n);
  for (int i = 0; i < m; ++i)
  {
    sol._F.setRowLabel(i, _F_point_unclustered.getRowLabel(i));
    for (int j = 0; j < n; ++j)
    {
      double f_ij = _cplex.getValue(_f[i][mutationSolVec[j]], solIdx);
      if (!g_tol.nonZero(f_ij))
      {
        f_ij = 0;
      }
      sol._F.set(i, j, f_ij);
    }
  }
  
  char buf[1024];
  for (int j = 0; j < n; ++j)
  {
    const StlIntVector& M = _toUnclusteredColumn[mutationSolVec[j]];
    std::string label;
    bool first = true;
    for (StlIntVectorIt it = M.begin(); it != M.end(); ++it)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        label += ";";
      }
      snprintf(buf, 1024, "%d", *it);
      label += buf;
    }
    
    sol._F.setColLabel(j, label);
  }
  
  sol._U = sol._T.getU(sol._F);
}
  
void IntMaxIlpSolver::initVariables()
{
  const int m = _F_interval_clustered.getNrRows();
  const int n = _F_interval_clustered.getNrCols();
  
  assert(_nodeCount == n);
  
  const Digraph& G = _G.getG();
  
#ifdef DEBUG
  char buf[1024];
#endif
  
  _x = IloBoolVarArray(_env, _arcCount + _nodeCount);
  int i = 0;
  for (NodeIt v(G); v != lemon::INVALID; ++v, ++i)
  {
#ifdef DEBUG
    snprintf(buf, 1024, "x_r_%d", i);
    _x[i] = IloBoolVar(_env, buf);
#else
    _x[i] = IloBoolVar(_env);
#endif
    _nodeToRootArcIndex[v] = i;
  }
  
  for (ArcIt a(G); a != lemon::INVALID; ++a, ++i)
  {
    _arcToIndex[a] = i;
    _indexToArc[i] = a;
    
#ifdef DEBUG
    snprintf(buf, 1024, "x_%d_%d",
             _nodeToIndex[G.source(a)],
             _nodeToIndex[G.target(a)]);
    _x[i] = IloBoolVar(_env, buf);
#else
    _x[i] = IloBoolVar(_env);
#endif
  }
  
  _f = IloNumVarMatrix(_env, m);
  _fx = IloNumVar3Matrix(_env, m);
  for (int i = 0; i < m; ++i)
  {
    _f[i] = IloNumVarArray(_env, n, 0, 0.5);
    _fx[i] = IloNumVarMatrix(_env, n);
    for (int j = 0; j < n; ++j)
    {
      const RealInterval& interval = _F_interval_clustered(i, j);
#ifdef DEBUG
      snprintf(buf, 1024, "f_%d_%d", i, j);
      _f[i][j].setName(buf);
#endif
      _f[i][j].setLB(interval.first);
      _f[i][j].setUB(interval.second);
      
      _fx[i][j] = IloNumVarArray(_env, _nodeCount + _arcCount, 0, 0.5);
#ifdef DEBUG
      for (int k = 0; k < _nodeCount + _arcCount; ++k)
      {
        Arc a = _indexToArc[k];
        if (a != lemon::INVALID)
        {
          Node u = G.source(a);
          Node v = G.target(a);

          snprintf(buf, 1024, "fx_%d_%d_(%d_%d)",
                   i, j,
                   _nodeToIndex[u], _nodeToIndex[v]);
          _fx[i][j][k].setName(buf);
        }
        else
        {
          snprintf(buf, 1024, "fx_%d_%d_(r_%d)",
                   i, j, k);
          _fx[i][j][k].setName(buf);
        }
      }
#endif
    }
  }
  
  _g = IloNumVarMatrix(_env, m);
  int org_n = _F_point_unclustered.getNrCols();
  for (int i = 0; i < m; ++i)
  {
    _g[i] = IloNumVarArray(_env, org_n, 0, 0.5);
#ifdef DEBUG
    for (int j = 0; j < org_n; ++j)
    {
      snprintf(buf, 1024, "g_%d_%d", i, j);
      _g[i][j].setName(buf);
    }
#endif
  }
}
  
void IntMaxIlpSolver::initConstraints()
{
  const Digraph& G = _G.getG();
  const int m = _F_interval_clustered.getNrRows();
  const int n = _F_interval_clustered.getNrCols();
  
  IloExpr sum(_env);
  
  // arborescence
  for (ArcIt jk(G); jk != lemon::INVALID; ++jk)
  {
    Node j = G.source(jk);
    sum += _x[_nodeToRootArcIndex[j]];
    for (InArcIt ij(G, j); ij != lemon::INVALID; ++ij)
    {
      sum += _x[_arcToIndex[ij]];
    }
    
    _model.add(sum >= _x[_arcToIndex[jk]]);
    sum.clear();
  }
  
  // at most one incoming arc per node
  for (NodeIt v(G); v != lemon::INVALID; ++v)
  {
    for (InArcIt a(G, v); a != lemon::INVALID; ++a)
    {
      sum += _x[_arcToIndex[a]];
    }
    _model.add(sum <= 1);
    sum.clear();
  }
  
  // there is exactly one root
  for (NodeIt v(G); v != lemon::INVALID; ++v)
  {
    sum += _x[_nodeToRootArcIndex[v]];
  }
  _model.add(sum == 1);
  sum.clear();
  
  // product
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      Node v_j = _indexToNode[j];
      for (int kl = 0; kl < _nodeCount + _arcCount; ++kl)
      {
        if (kl < _nodeCount)
        {
          // fake arc
          if (kl != _nodeToRootArcIndex[v_j])
          {
            _model.add(_fx[i][j][kl] == 0);
          }
          else
          {
            _model.add(_fx[i][j][kl] <= _f[i][j]);
            _model.add(_fx[i][j][kl] <= _x[kl]);
            _model.add(_fx[i][j][kl] >= _f[i][j] + _x[kl] - 1);
          }
        }
        else
        {
          Arc a = _indexToArc[kl];
          assert(a != lemon::INVALID);
          Node v_k = G.source(a);
          Node v_l = G.target(a);
          
          if (v_k == v_j || v_l == v_j)
          {
            _model.add(_fx[i][j][kl] <= _f[i][j]);
            _model.add(_fx[i][j][kl] <= _x[kl]);
            _model.add(_fx[i][j][kl] >= _f[i][j] + _x[kl] - 1);
          }
          else
          {
            _model.add(_fx[i][j][kl] == 0);
          }
        }
      }
    }
  }
  
  // sum rule constraint
  for (NodeIt v_k(G); v_k != lemon::INVALID; ++v_k)
  {
    int k = _nodeToIndex[v_k];
    for (int i = 0; i < m; ++i)
    {
      sum += _fx[i][k][_nodeToRootArcIndex[v_k]];
      for (InArcIt a(G, v_k); a != lemon::INVALID; ++a)
      {
//        Node v_j = G.source(a);
//        int j = _nodeToIndex[v_j];
        sum += _fx[i][k][_arcToIndex[a]];
      }
      for (OutArcIt a(G, v_k); a != lemon::INVALID; ++a)
      {
        Node v_l = G.target(a);
        int l = _nodeToIndex[v_l];
        sum -= _fx[i][l][_arcToIndex[a]];
      }
      
      _model.add(sum >= 0);
      sum.clear();
    }
  }
  
  // VAF deviation
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      const StlIntVector& M = _toUnclusteredColumn[j];
      for (StlIntVectorIt it = M.begin(); it != M.end(); ++it)
      {
        // we need to truncate the point estimate f at 0.5 because g <= 0.5
        double f = std::min(0.5, _F_point_unclustered(i, *it));
        _model.add(_g[i][*it] >= f - _f[i][j]);
        _model.add(_g[i][*it] >= _f[i][j] - f);
      }
    }
  }
  
  sum.end();
}
  
void IntMaxIlpSolver::initObjective()
{
  IloExpr sum(_env);
  
  for (int i = 0; i < _nodeCount + _arcCount; ++i)
  {
    sum += _x[i];
  }
  
  const int m = _F_point_unclustered.getNrRows();
  const int org_n = _F_point_unclustered.getNrCols();
  const double frac = 1.0 / (m * org_n);
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < org_n; ++j)
    {
      sum -= frac * _g[i][j];
    }
  }
  
  _model.add(IloObjective(_env, sum, IloObjective::Maximize));
  
  sum.end();
}
  
}