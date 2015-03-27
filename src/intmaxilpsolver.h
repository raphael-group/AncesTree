/*
 *  intmaxilpsolver.h
 *
 *   Created on: 8-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef INTMAXILPSOLVER_H
#define INTMAXILPSOLVER_H

#include <ilcplex/ilocplex.h>
#include "baseancestrygraph.h"
#include "realintervalmatrix.h"
#include "maxsolution.h"
#include <vector>

namespace vaff {

class IntMaxIlpSolver
{
public:
  IntMaxIlpSolver(const BaseAncestryGraph& G,
                  const RealIntervalMatrix& F_interval_clustered,
                  const RealMatrix& F_point_unclustered,
                  const StlIntMatrix& toUnclusteredColumn,
                  int timeLimit);
  
  ~IntMaxIlpSolver();
  
  bool solve(MaxSolution& solution);
  
  void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
  }
  
protected:
  typedef BaseAncestryGraph::Digraph Digraph;
  DIGRAPH_TYPEDEFS(Digraph);
  
  typedef std::vector<Arc> ArcVector;
  typedef std::vector<Node> NodeVector;
  
  void initVariables();
  void initConstraints();
  void initObjective();
  
  void printVariables(int solIdx, std::ostream& out) const;
  
  void processSolution(int solIdx, MaxSolution::Triple& sol) const;
  
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolArray> IloBoolMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumArray> IloNumMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  
  const BaseAncestryGraph& _G;
  const RealIntervalMatrix& _F_interval_clustered;
  const RealMatrix& _F_point_unclustered;
  const StlIntMatrix& _toUnclusteredColumn;
  const int _timeLimit;
  
  const int _nodeCount;
  const IntNodeMap& _nodeToIndex;
  const NodeVector& _indexToNode;
  IntNodeMap _nodeToRootArcIndex;
  
  const int _arcCount;
  IntArcMap _arcToIndex;
  ArcVector _indexToArc;
  
  IloEnv _env;
  IloModel _model;
  IloCplex _cplex;
  
  // x[(j,k)] : arc (v_j,v_k) is in the solution
  IloBoolVarArray _x;
  // f[i][j] : corrected VAF for sample i and mutation j
  IloNumVarMatrix _f;
  // fx[i][j][(k,l)] = f[i][j] * x[(k,l)]
  IloNumVar3Matrix _fx;
  // f[i][j] : VAF deviation for sample i and unclustered mutation j
  IloNumVarMatrix _g;
};
  
}

#endif // MAXILPSOLVER_H