/*
 * =====================================================================================
 *
 *       Filename:  polyfit.hpp
 *
 *    Description:  Source code to handle polynomial fits using least squares and to
 *                  find extrema of quadratic and cubic polynomials.
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

#ifndef POLYFIT_HPP
#define POLYFIT_HPP

#include "matrix.hpp"
#include "verbosity.hpp"

/** Returns a vector of coefficients of a polynomial of degree n that is fit to
 * the data.
 */
vector_t polyfit(const vector_t &x, const vector_t &y, size_t n, Verbosity verbosity);

double quadratic_extremum(const vector_t &x, const vector_t &y, Verbosity verbosity);
std::pair<double,double> cubic_extrema(const vector_t &x, const vector_t &y, Verbosity verbosity);

/** Find the maximum or minimum of the quadratic polynomial c[2]*x^2 + c[1]*x + c[0] = 0 */
double quadratic_extremum(const vector_t &c, Verbosity verbose);

/** Find the minimum and maximum of the cubic polynomial c[2]*x^2 + c[1]*x + c[0] = 0
 * Return: a pair of values with the first one being the minimum, and the second being the maximum.
 */
std::pair<double,double> cubic_extrema(const vector_t &c, Verbosity verbose);

/** Find the roots of the quadratic polynomial a*x^2 + b*x + c = 0
 * Return: a pair of value for the two roots, with the smaller one in the first element.
 * Warning: No checking is done for e.g. single or no roots.
 */
std::pair<double,double> solve_quadratic(double a, double b, double c, Verbosity verbose);

/** Find the roots of the quadratic polynomial c[2]*x^2 + c[1]*x + c[0] = 0
 * Return: a pair of value for the two roots, with the smaller one in the first element.
 * Warning: No checking is done for e.g. single or no roots.
 */
std::pair<double,double> solve_quadratic(const vector_t &c, Verbosity verbose);

vector_t interpolate3(double x1, double x2, double f1, double f2, double g1, double g2);
vector_t interpolate2(double x1, double x2, double f1, double f2, double g1, Verbosity verbose);
vector_t interpolate2(double x1, double x2, double g1, double g2);

#endif

