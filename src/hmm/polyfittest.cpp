/*
 * =====================================================================================
 *
 *       Filename:  polyfittest.cpp
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

using namespace std;

int main(int argc, const char **argv)
{
  size_t k = 2;
  if(argc == 2)
    k = atoi(argv[1]);
  size_t n = 10;
  vector_t x(n);
  vector_t y(n);
  for(size_t i = 0; i < n; i++)
    x(i) = i;
  for(size_t i = 0; i < n; i++)
    y(i) = i*i;

  vector_t coeff = polyfit(x, y, k, Verbosity::info);
  double mq = quadratic_extremum(x, y, Verbosity::info);
  std::pair<double,double> mc = cubic_extrema(x, y, Verbosity::info);

  std::cout << "x =";
  for(auto &u: x)
    std::cout << " " << u;
  std::cout << std::endl;

  std::cout << "y =";
  for(auto &u: y)
    std::cout << " " << u;
  std::cout << std::endl;

  std::cout << "alpha = ";
  for(auto &u: coeff)
    std::cout << " " << u;
  std::cout << std::endl;

  std::cout << "The inferred minimum of the quadratic is at " << mq << std::endl;
  std::cout << "The inferred minimum of the cubic is at " << mc.first << std::endl;
  std::cout << "The inferred maximum of the cubic is at " << mc.second << std::endl;
  
  return EXIT_SUCCESS;
}

