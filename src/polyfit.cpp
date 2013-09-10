/*
 * =====================================================================================
 *
 *       Filename:  polyfit.cpp
 *
 *    Description:  Source code to compute fits of polynomial using least squares
 *
 *        Version:  1.0
 *        Created:  09/27/2012 03:32:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include "polyfit.hpp"
#include "matrix_inverse.hpp"

double quadratic_extremum(const vector_t &coeff, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    std::cout << "Quadratic extremum of " << coeff;
  double x = -coeff(1) / (2*coeff(2));
  if(verbosity >= Verbosity::debug) {
    std::cout << " = " << x << std::endl;
    std::cout << "It's a ";
    if(coeff(2) > 0)
      std::cout << "minimum." << std::endl;
    else if(coeff(2) < 0)
      std::cout << "maximum." << std::endl;
    else
      std::cout << "flat line!" << std::endl;
  }
  return(x);
}

std::pair<double,double> solve_quadratic(double a, double b, double c, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    std::cout << "Solve quadratic for a = " << a << " b = " << b << " c = " << c;
  std::pair<double,double> x;
  x.first = (-b - sqrt(b*b-4*a*c)) / (2*a);
  x.second = (-b + sqrt(b*b-4*a*c)) / (2*a);
  if(x.first > x.second) { // make sure the first element is the smaller one
    if(verbosity >= Verbosity::debug)
      std::cout << " (swapping!)";
    std::swap(x.first, x.second);
  }
  if(verbosity >= Verbosity::debug)
    std::cout << " -> " << x.first << " and " << x.second << std::endl;
  return(x);
}

std::pair<double,double> solve_quadratic(const vector_t &coeff, Verbosity verbosity) {
  return(solve_quadratic(coeff(2), coeff(1), coeff(0), verbosity));
}

std::pair<double,double> cubic_extrema(const vector_t &coeff, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    std::cout << "Cubic extrema for " << coeff << std::endl;
  assert(coeff.size() == 4);
  assert(coeff[3] != 0);
  vector_t deriv(3);
  for(size_t i = 0; i < 3; i++)
    deriv(i) = (i+1) * coeff(i+1);

  std::pair<double, double> extrema = solve_quadratic(deriv, verbosity);

  if(coeff(3) > 0)
    std::swap(extrema.first, extrema.second);

  if(verbosity >= Verbosity::debug)
    std::cout << "Cubic extrema at " << extrema.first << " " << extrema.second << std::endl;
  return(extrema);
}

double quadratic_extremum(const vector_t &x, const vector_t &y, Verbosity verbosity) {
  return(quadratic_extremum(polyfit(x, y, 2, verbosity), verbosity)); 
}

std::pair<double,double> cubic_extrema(const vector_t &x, const vector_t &y, Verbosity verbosity) {
  vector_t coeff = polyfit(x, y, 3, verbosity); 
  return(cubic_extrema(coeff, verbosity));
}

/** Returns a vector of coefficients of a polynomial of degree n that is fit to the data */
vector_t polyfit(const vector_t &x, const vector_t &y, size_t n, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    std::cout << "Fitting " << n << "-th order polynomial." << std::endl;
  if(n >= x.size()) {
    std::cout << "Error in polynomial fitting: cannot fit " << n << "-th order polynomial to " << x.size() << " data points." << std::endl;
    exit(-1);
  }

  using namespace boost::numeric::ublas;
  matrix_t X(x.size(), n+1);
  for(size_t i = 0; i < x.size(); i++)
    for(size_t k = 0; k < n+1; k++)
      X(i,k) = pow(x(i), k);

  if(verbosity >= Verbosity::debug)
    std::cout << "X = " << X << std::endl;

  matrix_t P = prod(trans(X), X);

  if(verbosity >= Verbosity::debug)
    std::cout << "P = " << P << std::endl;

  matrix_t Q(n+1,n+1);
  bool inverted = InvertMatrix(P, Q);

  if(not inverted) {
    std::cout << "Couldn't invert matrix." << std::endl;
    exit(-1);
  }
  if(verbosity >= Verbosity::debug)
    std::cout << "Q = " << Q << std::endl;

  vector_t coeff = prod(prod(Q, trans(X)), y);
  return(coeff);
}

vector_t interpolate3(double x1, double x2, double f1, double f2, double g1, double g2) {
  matrix_t P(4, 4);
  P(0,0) = 1;
  P(1,0) = 1;
  P(2,0) = 0;
  P(3,0) = 0;
  P(0,1) = x1;
  P(1,1) = x2;
  P(2,1) = 1;
  P(3,1) = 1;
  P(0,2) = x1*x1;
  P(1,2) = x2*x2;
  P(2,2) = 2*x1;
  P(3,2) = 2*x2;
  P(0,3) = x1*x1*x1;
  P(1,3) = x2*x2*x2;
  P(2,3) = 3*x1*x1;
  P(3,3) = 3*x2*x2;
  matrix_t Q(4, 4);
  InvertMatrix(P, Q);
  vector_t A(4);
  A(0) = f1;
  A(1) = f2;
  A(2) = g1;
  A(3) = g2;
  vector_t coeff = boost::numeric::ublas::prod(Q, A);
//  std::cout << P << std::endl << Q << std::endl << A << std::endl;
//  std::cout << coeff << std::endl;
  return(coeff);
}

vector_t interpolate2(double x1, double x2, double f1, double f2, double g1, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    std::cout << "interpolate2 " << x1 << " " << x2 << " " << f1 << " " << f2 << " " << g1 << std::endl;
  double a = - (f1 - f2 - g1 * (x1 - x2)) / ((x1 - x2) * (x1 - x2));
  double b = g1 - 2 * a * x1;
  double c = f1 - a * x1 * x1 - b * x1;
  vector_t coeff(3);
  coeff(0) = c;
  coeff(1) = b;
  coeff(2) = a;
  return(coeff);
}

vector_t interpolate2(double x1, double x2, double g1, double g2) {
  vector_t coeff(3);
  coeff(2) = 0.5 * (g1 - g2) / (x1 - x2);
  coeff(1) = g1 - 2 * coeff(2) * x1;
  coeff(0) = 0;
  return(coeff);
}
