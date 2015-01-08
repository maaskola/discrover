/*
 * =====================================================================================
 *
 *       Filename:  polyfittest.cpp
 *
 *    Description:  Compute fits of polynomial using least squares
 *
 *        Created:  09/27/2012 03:32:29 PM
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */
#include <stdlib.h>

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include "polyfit.hpp"

using namespace std;

int main(int argc, const char **argv) {
  size_t k = 2;
  if (argc == 2)
    k = atoi(argv[1]);
  size_t n = 10;
  vector_t x(n);
  vector_t y(n);
  for (size_t i = 0; i < n; i++)
    x(i) = i;
  for (size_t i = 0; i < n; i++)
    y(i) = i * i;

  vector_t coeff = polyfit(x, y, k, Verbosity::info);
  double mq = quadratic_extremum(x, y, Verbosity::info);
  pair<double, double> mc = cubic_extrema(x, y, Verbosity::info);

  cout << "x =";
  for (auto &u : x)
    cout << " " << u;
  cout << endl;

  cout << "y =";
  for (auto &u : y)
    cout << " " << u;
  cout << endl;

  cout << "alpha = ";
  for (auto &u : coeff)
    cout << " " << u;
  cout << endl;

  cout << "The inferred minimum of the quadratic is at " << mq << endl;
  cout << "The inferred minimum of the cubic is at " << mc.first << endl;
  cout << "The inferred maximum of the cubic is at " << mc.second << endl;

  return EXIT_SUCCESS;
}
