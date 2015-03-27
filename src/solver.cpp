/*
 *  solver.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "solver.h"

namespace vaff {

Solver::Solver(const RealMatrix& F)
  : _F(F)
{
}

Solver::~Solver()
{
}

} // namespace vaff
