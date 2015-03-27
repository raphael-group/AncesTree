/*
 *  comparison.h
 *
 *   Created on: 12-mar-2015
 *       Author: M. El-Kebir
 */

#ifndef COMPARISON_H
#define COMPARISON_H

#include "utils.h"
#include "maxsolution.h"
#include "clonaltree.h"

#include <set>
#include <algorithm>

namespace vaff {
  
class Comparison
{
public:
  Comparison(const MaxSolution& referenceSolution,
             const MaxSolution::Triple& solution,
             const MaxSolution::Triple& whitelistSolution);
  
  double coverage() const;
  
  double deltaF() const;
  
  double deltaU() const;
  
  void accuracyB(double& clustered,
                 double& ancestral,
                 double& incomparable) const;
  
  void recallB(double& clustered,
               double& ancestral,
               double& incomparable) const;

private:
  typedef std::pair<int, int> StlIntPair;
  typedef std::set<StlIntPair> StlIntSet;
  
  const MaxSolution& _referenceSolution;
  const MaxSolution::Triple& _solution;
  StlIntMatrix _toMutationsFromSol;
  StlIntMatrix _toMutationsFromRef;
  StlIntVector _toSolCluster;
  StlIntVector _toRefCluster;
  StlBoolVector _whitelist;
  
  void constructMappings(const MaxSolution::Triple& whitelistSolution);
  
  void determinePairs(const MaxSolution::Triple& solution,
                      const StlIntMatrix& toMutations,
                      const StlIntVector& toCluster,
                      StlIntSet& P_clustered,
                      StlIntSet& N_clustered,
                      StlIntSet& P_ancestral,
                      StlIntSet& N_ancestral,
                      StlIntSet& P_incomparable,
                      StlIntSet& N_incomparable) const;
  
  void determineNodeMapping(StlIntVector& ref2sol,
                            StlIntVector& sol2ref) const;
  
  double recall(const StlIntSet& S_ref,
                const StlIntSet& S_sol) const;
  
  double accuracy(const StlIntSet& P_ref,
                  const StlIntSet& N_ref,
                  const StlIntSet& P_sol,
                  const StlIntSet& N_sol) const;
};
  
} // namespace vaff

#endif // COMPARISON_H
