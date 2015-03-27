/*
 *  probancestrygraph.h
 *
 *   Created on: 10-jan-2015
 *       Author: M. El-Kebir
 */

#ifndef PROBANCESTRYGRAPH_H
#define PROBANCESTRYGRAPH_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <lemon/connectivity.h>
#include "baseancestrygraph.h"
#include "utils.h"
#include "ancestrymatrix.h"
#include "readcountmatrix.h"

namespace vaff {
  
class ProbAncestryGraph : public BaseAncestryGraph
{
public:
  typedef std::vector<std::string> StringVector;
  typedef StringVector::const_iterator StringVectorIt;
  
  typedef lemon::ListDigraph Digraph;
  DIGRAPH_TYPEDEFS(Digraph);
  typedef std::vector<Node> NodeVector;
  
  ProbAncestryGraph();
  
  ProbAncestryGraph(const AncestryMatrix& A,
                    const ReadCountMatrix& R,
                    double alpha,
                    double gamma);
  
  int numberOfNodesInfCI() const;
  
  virtual void removeCycles(const AncestryMatrix& A,
                            double alpha,
                            StlIntMatrix& toOriginalColumns) const;
  
  void removeCyclesYoshiko(const AncestryMatrix& A,
                           double alpha,
                           StlIntMatrix& toOriginalColumns) const;
  
  void contract(const AncestryMatrix& A,
                const StlIntMatrix& toOrginalColumns,
                double beta,
                ProbAncestryGraph& H) const;
  
  void writeDOT(const ReadCountMatrix& R,
                std::ostream& out) const;
  
  void writeDOT(const ReadCountMatrix& R,
                const StlIntMatrix& M,
                std::ostream& out) const;
  
  void writeDOT(const ReadCountMatrix& R,
                const StlIntVector& S,
                std::ostream& out) const;
  
  static void createLabels(const StlIntMatrix& toOrginalColumns,
                           StringVector& label);
  
  static void printClusterSize(const StlIntMatrix& toOrginalColumns,
                               std::ostream& out);
  
  static void printInterClusterCoherence(const Digraph& G,
                                         const IntNodeMap& nodeToMutation,
                                         const StringVector& mutationToLabel,
                                         const AncestryMatrix& A,
                                         std::ostream& out);
  
  static double interClusterCoherence(const Digraph& G,
                                      const IntNodeMap& nodeToMutation,
                                      const StringVector& mutationToLabel,
                                      const AncestryMatrix& A);
  
  static void printIntraClusterCoherence(const Digraph& G,
                                         const IntNodeMap& nodeToMutation,
                                         const StringVector& mutationToLabel,
                                         const AncestryMatrix& A,
                                         std::ostream& out);
  
  static void intraClusterCoherence(const Digraph& G,
                                    const IntNodeMap& nodeToMutation,
                                    const StringVector& mutationToLabel,
                                    const AncestryMatrix& A,
                                    double& incomparable,
                                    double& ancestral,
                                    double& comparable);
  
protected:
  BoolArcMap _intermediateArc;
  
private:
  void constructGraph(const AncestryMatrix& A,
                      const ReadCountMatrix& R,
                      double alpha,
                      double gamma);
};
  
} // namespace vaff

#endif // PROBANCESTRYGRAPH_H
