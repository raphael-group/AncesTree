/*
 *  utils.h
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <string>
#include <istream>
#include <lemon/tolerance.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

namespace vaff {

typedef std::vector<bool> StlBoolVector;
typedef std::vector<StlBoolVector> StlBoolMatrix;
  
typedef std::vector<int> StlIntVector;
typedef StlIntVector::const_iterator StlIntVectorIt;
typedef std::vector<StlIntVector> StlIntMatrix;
typedef StlIntMatrix::const_iterator StlIntMatrixIt;
  
typedef std::vector<double> StlDoubleVector;
typedef std::vector<StlDoubleVector> StlDoubleMatrix;
  
typedef std::pair<double, double> RealInterval;
typedef std::vector<RealInterval> StlRealIntervalVector;
typedef std::vector<StlRealIntervalVector> StlRealIntervalMatrix;
  
typedef std::pair<int, int> IntPair;
  
std::ostream& operator<<(std::ostream& out, const StlBoolMatrix& M);
std::istream& operator>>(std::istream& in, StlBoolMatrix& M);

std::ostream& operator<<(std::ostream& out, const StlDoubleMatrix& M);
std::istream& operator>>(std::istream& in, StlDoubleMatrix& M);
  
std::ostream& operator<<(std::ostream& out, const StlRealIntervalMatrix& M);
std::istream& operator>>(std::istream& in, StlRealIntervalMatrix& M);

StlBoolVector discretize(const StlDoubleVector& v);

StlBoolMatrix discretize(const StlDoubleMatrix& M);
  
std::istream& getline(std::istream& is, std::string& t);
  
extern lemon::Tolerance<double> g_tol;
  
} // namespace vaff

#endif // UTILS_H
