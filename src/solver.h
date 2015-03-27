/*
 *  solver.h
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "utils.h"
#include "solution.h"

namespace vaff {

class Solver
{
public:
  Solver(const RealMatrix& V);
  
  virtual ~Solver();
  
  virtual bool solve(Solution& solution) = 0;
  
protected:
  const RealMatrix& _F;
};
  
} // namespace vaff

#endif // SOLVER_H




