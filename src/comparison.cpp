/*
 *  comparison.cpp
 *
 *   Created on: 12-mar-2015
 *       Author: M. El-Kebir
 */

#include "comparison.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <lemon/matching.h>

namespace vaff {

Comparison::Comparison(const MaxSolution& referenceSolution,
                       const MaxSolution::Triple& solution,
                       const MaxSolution::Triple& whitelistSolution)
  : _referenceSolution(referenceSolution)
  , _solution(solution)
  , _toMutationsFromSol()
  , _toMutationsFromRef()
  , _toSolCluster()
  , _toRefCluster()
{
  constructMappings(whitelistSolution);
}
  
void Comparison::constructMappings(const MaxSolution::Triple& whitelistSolution)
{
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;

  int n = _solution._F.getNrCols();
  
  _toMutationsFromSol = StlIntMatrix(n, StlIntVector());
  _toSolCluster = StlIntVector(_referenceSolution.getF().getNrCols(), -1);
  
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, _solution._F.getColLabel(j), boost::is_any_of(";"));
    
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      if (*it2 != "")
      {
        int org_j = boost::lexical_cast<int>(*it2);
        
        _toMutationsFromSol[j].push_back(org_j);
        _toSolCluster[org_j] = j;
      }
    }
  }
  
  const RealMatrix& refF = _referenceSolution.solution(0)._F;
  n = _referenceSolution.solution(0)._F.getNrCols();
  
  _toMutationsFromRef = StlIntMatrix(n, StlIntVector());
  _toRefCluster = StlIntVector(_referenceSolution.getF().getNrCols(), -1);
  
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, refF.getColLabel(j), boost::is_any_of(";"));
    
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      int org_j = boost::lexical_cast<int>(*it2);
      
      _toMutationsFromRef[j].push_back(org_j);
      _toRefCluster[org_j] = j;
    }
  }
  
  n = whitelistSolution._F.getNrCols();
  _whitelist = StlBoolVector(_referenceSolution.getF().getNrCols(), false);
  for (int j = 0; j < n; ++j)
  {
    StringVector s;
    boost::split(s, whitelistSolution._F.getColLabel(j), boost::is_any_of(";"));
    
    for (StringVectorIt it2 = s.begin(); it2 != s.end(); ++it2)
    {
      if (*it2 != "")
      {
        int org_j = boost::lexical_cast<int>(*it2);
        _whitelist[org_j] = true;
      }
    }
  }
}
  
double Comparison::coverage() const
{
  int sol_mut_count = 0;
  for (StlIntVectorIt sol_it = _toSolCluster.begin();
       sol_it != _toSolCluster.end(); ++sol_it)
  {
    if (*sol_it != -1)
      ++sol_mut_count;
  }
  
  int ref_mut_count = 0;
  for (StlIntVectorIt ref_it = _toRefCluster.begin();
       ref_it != _toRefCluster.end(); ++ref_it)
  {
    if (*ref_it != -1)
      ++ref_mut_count;
  }
  
  return (double) sol_mut_count / (double) ref_mut_count;
}
  
double Comparison::deltaF() const
{
  const RealMatrix& F = _solution._F;
  const RealMatrix& refF = _referenceSolution.getF();
  
  int n = F.getNrCols();
  int m = F.getNrRows();
  
  assert(m == refF.getNrRows());
  
  double delta = 0;
  int entries = 0;
  for (int j = 0; j < n; ++j)
  {
    const StlIntVector& M = _toMutationsFromSol[j];
    for (StlIntVectorIt it = M.begin(); it != M.end(); ++it)
    {
      int org_j = *it;
      for (int p = 0; p < m; ++p)
      {
//        std::cout << p << "\t" << org_j << ":\t|"
//                  << F(p, j) << " - " << refF(p, org_j)
//                  << "| = " << fabs(F(p, j) - refF(p, org_j)) << std::endl;
        delta += fabs(F(p, j) - refF(p, org_j));
        ++entries;
      }
    }
  }
  
  return delta / entries;
}
  
double Comparison::deltaU() const
{
  double delta = 0;
  int entries = 0;
  
  StlIntVector ref2sol;
  StlIntVector sol2ref;
  
  determineNodeMapping(ref2sol, sol2ref);
  
  const RealMatrix& solU = _solution._U;
  const RealMatrix& refU = _referenceSolution.solution(0)._U;
  
  assert(solU.getNrRows() == refU.getNrRows());
  
  int m = refU.getNrRows();
  int n = refU.getNrCols();
  
  for (int j = 0; j < n; ++j)
  {
    if (ref2sol[j] != -1)
    {
      for (int p = 0; p < m; ++p)
      {
//        std::cout << p << "\t" << j << "/" << ref2sol[j] << ":\t|"
//                  << refU(p, j) << " - " << solU(p, ref2sol[j])
//                  << "| = " << fabs(refU(p, j) - solU(p, ref2sol[j])) << std::endl;
        
        delta += fabs(refU(p, j) - solU(p, ref2sol[j]));
        ++entries;
      }
    }
  }
  
  return delta / entries;
}
  
void Comparison::recallB(double& clustered,
                         double& ancestral,
                         double& incomparable) const
{
  StlIntSet P_ref_clustered, N_ref_clustered;
  StlIntSet P_ref_ancestral, N_ref_ancestral;
  StlIntSet P_ref_incomparable, N_ref_incomparable;
  
  determinePairs(_referenceSolution.solution(0),
                 _toMutationsFromRef,
                 _toRefCluster,
                 P_ref_clustered, N_ref_clustered,
                 P_ref_ancestral, N_ref_ancestral,
                 P_ref_incomparable, N_ref_incomparable);
  
  StlIntSet P_sol_clustered, N_sol_clustered;
  StlIntSet P_sol_ancestral, N_sol_ancestral;
  StlIntSet P_sol_incomparable, N_sol_incomparable;
  
  determinePairs(_solution,
                 _toMutationsFromSol,
                 _toSolCluster,
                 P_sol_clustered, N_sol_clustered,
                 P_sol_ancestral, N_sol_ancestral,
                 P_sol_incomparable, N_sol_incomparable);
  
  clustered = recall(P_ref_clustered, P_sol_clustered);
  ancestral = recall(P_ref_ancestral, P_sol_ancestral);
  incomparable = recall(P_ref_incomparable, P_sol_incomparable);
}
  
double Comparison::recall(const StlIntSet& S_ref,
                          const StlIntSet& S_sol) const
{
  StlIntSet intersection;
  std::set_intersection(S_ref.begin(), S_ref.end(),
                        S_sol.begin(), S_sol.end(),
                        std::inserter(intersection, intersection.begin()));
  
  return (double) intersection.size() / (double) S_ref.size();
}
  
void Comparison::accuracyB(double& clustered,
                           double& ancestral,
                           double& incomparable) const
{
  StlIntSet P_ref_clustered, N_ref_clustered;
  StlIntSet P_ref_ancestral, N_ref_ancestral;
  StlIntSet P_ref_incomparable, N_ref_incomparable;
  
  determinePairs(_referenceSolution.solution(0),
                 _toMutationsFromRef,
                 _toRefCluster,
                 P_ref_clustered, N_ref_clustered,
                 P_ref_ancestral, N_ref_ancestral,
                 P_ref_incomparable, N_ref_incomparable);
  
  StlIntSet P_sol_clustered, N_sol_clustered;
  StlIntSet P_sol_ancestral, N_sol_ancestral;
  StlIntSet P_sol_incomparable, N_sol_incomparable;
  
  determinePairs(_solution,
                 _toMutationsFromSol,
                 _toSolCluster,
                 P_sol_clustered, N_sol_clustered,
                 P_sol_ancestral, N_sol_ancestral,
                 P_sol_incomparable, N_sol_incomparable);
  
  clustered = accuracy(P_ref_clustered, N_ref_clustered,
                       P_sol_clustered, N_sol_clustered);
  ancestral = accuracy(P_ref_ancestral, N_ref_ancestral,
                       P_sol_ancestral, N_sol_ancestral);
  incomparable = accuracy(P_ref_incomparable, N_ref_incomparable,
                          P_sol_incomparable, N_sol_incomparable);
}
  
double Comparison::accuracy(const StlIntSet& P_ref,
                            const StlIntSet& N_ref,
                            const StlIntSet& P_sol,
                            const StlIntSet& N_sol) const
{
  assert(P_ref.size() + N_ref.size() == P_sol.size() + N_sol.size());
  
  StlIntSet TP;
  std::set_intersection(P_ref.begin(), P_ref.end(),
                        P_sol.begin(), P_sol.end(),
                        std::inserter(TP, TP.begin()));
  
  StlIntSet TN;
  std::set_intersection(N_ref.begin(), N_ref.end(),
                        N_sol.begin(), N_sol.end(),
                        std::inserter(TN, TN.begin()));
  
  int total = P_ref.size() + N_ref.size();
  
  return (double) (TP.size() + TN.size()) / (double) total;
}
  
void Comparison::determinePairs(const MaxSolution::Triple& solution,
                                const StlIntMatrix& toMutations,
                                const StlIntVector& toCluster,
                                StlIntSet& P_clustered,
                                StlIntSet& N_clustered,
                                StlIntSet& P_ancestral,
                                StlIntSet& N_ancestral,
                                StlIntSet& P_incomparable,
                                StlIntSet& N_incomparable) const
{
  P_clustered.clear();
  N_clustered.clear();
  P_ancestral.clear();
  N_ancestral.clear();
  P_incomparable.clear();
  N_incomparable.clear();
  
  int org_j = 0;
  for (StlIntVectorIt it_1 = toCluster.begin();
       it_1 != toCluster.end(); ++it_1, ++org_j)
  {
    int j = *it_1;
    if (j == -1) continue;
    
    int org_k = org_j + 1;
    for (StlIntVectorIt it_2 = it_1 + 1;
         it_2 != toCluster.end(); ++it_2, ++org_k)
    {
      int k = *it_2;
      if (k == -1) continue;
      
      if (!_whitelist[org_j] || !_whitelist[org_k])
        continue;
      
      if (j == k)
      {
        P_clustered.insert(std::make_pair(std::min(org_j, org_k),
                                          std::max(org_j, org_k)));
        N_ancestral.insert(std::make_pair(org_j, org_k));
        N_ancestral.insert(std::make_pair(org_k, org_j));
        N_incomparable.insert(std::make_pair(std::min(org_j, org_k),
                                             std::max(org_j, org_k)));
      }
      else
      {
        N_clustered.insert(std::make_pair(std::min(org_j, org_k),
                                          std::max(org_j, org_k)));
        
        if (solution._T.ancestral(j, k))
        {
          P_ancestral.insert(std::make_pair(org_j, org_k));
          N_ancestral.insert(std::make_pair(org_k, org_j));
          N_incomparable.insert(std::make_pair(std::min(org_j, org_k),
                                               std::max(org_j, org_k)));
        }
        else if (solution._T.ancestral(k, j))
        {
          P_ancestral.insert(std::make_pair(org_k, org_j));
          N_ancestral.insert(std::make_pair(org_j, org_k));
          N_incomparable.insert(std::make_pair(std::min(org_j, org_k),
  std::max(org_j, org_k)));
        }
        else
        {
          assert(solution._T.incomparable(j, k));
          P_incomparable.insert(std::make_pair(std::min(org_j, org_k),
                                               std::max(org_j, org_k)));
          N_ancestral.insert(std::make_pair(org_j, org_k));
          N_ancestral.insert(std::make_pair(org_k, org_j));
        }
      }
    }
  }
}
  
void Comparison::determineNodeMapping(StlIntVector& ref2sol,
                                      StlIntVector& sol2ref) const
{
  const ClonalTree& refT = _referenceSolution.solution(0)._T;
  const ClonalTree& solT = _solution._T;
  
  int ref_n = refT.getB().getNrRows();
  int sol_n = solT.getB().getNrRows();
  int n = std::max(ref_n, sol_n);
  
  ref2sol = StlIntVector(ref_n, -1);
  sol2ref = StlIntVector(sol_n, -1);
  
  // make bipartite graph
  lemon::ListBpGraph G;
  lemon::ListBpGraph::RedNodeMap<int> redMut(G);
  lemon::ListBpGraph::BlueNodeMap<int> blueMut(G);
  lemon::ListBpGraph::EdgeMap<int> cost(G);
  
  // add nodes (and dummy nodes)
  for (int i = 0; i < n; ++i)
  {
    lemon::ListBpGraph::BlueNode b = G.addBlueNode();
    lemon::ListBpGraph::RedNode r = G.addRedNode();
    
    redMut[r] = i < ref_n ? i : -1;
    blueMut[b] = i < sol_n ? i : -1;
  }
  
  // add edges
  for (lemon::ListBpGraph::RedNodeIt r(G); r != lemon::INVALID; ++r)
  {
    int mut_r = redMut[r];
    for (lemon::ListBpGraph::BlueNodeIt b(G); b != lemon::INVALID; ++b)
    {
      int mut_b = blueMut[b];
      
      lemon::ListBpGraph::Edge e = G.addEdge(r, b);
      
      if (mut_r == -1 && mut_b == -1)
      {
        cost[e] = 0;
      }
      else if (mut_r == -1)
      {
        cost[e] = -_toMutationsFromSol[mut_b].size();
      }
      else if (mut_b == -1)
      {
        cost[e] = -_toMutationsFromRef[mut_r].size();
      }
      else
      {
        StlIntVector sym_diff;
        const StlIntVector& S_r = _toMutationsFromRef[mut_r];
        const StlIntVector& S_b = _toMutationsFromSol[mut_b];
        std::set_symmetric_difference(S_r.begin(), S_r.end(),
                                      S_b.begin(), S_b.end(),
                                      std::inserter(sym_diff, sym_diff.begin()));
        
        cost[e] = -sym_diff.size();
      }
    }
  }
  
  lemon::MaxWeightedPerfectMatching<lemon::ListBpGraph> mwpm(G, cost);
  mwpm.run();
  
//  std::cout << mwpm.matchingWeight() << std::endl;
  
  for (lemon::ListBpGraph::EdgeIt e(G); e != lemon::INVALID; ++e)
  {
    lemon::ListBpGraph::RedNode r = G.redNode(e);
    lemon::ListBpGraph::BlueNode b = G.blueNode(e);
    
    int mut_r = redMut[r];
    int mut_b = blueMut[b];
    
    if (mwpm.matching(e))
    {
      if (mut_r != -1)
      {
        ref2sol[mut_r] = mut_b;
      }
      if (mut_b != -1)
      {
        sol2ref[mut_b] = mut_r;
      }
//      std::cout << "*";
    }
    
//    std::cout << redMut[r] << " -- " << blueMut[b] << " : " << cost[e] << std::endl;
  }
}
  
} // namespace vaff
