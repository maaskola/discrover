/* =====================================================================================
 * Copyright (c) 2011, Jonas Maaskola
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 *
 *       Filename:  aux.hpp
 *
 *    Description:  Header for auxiliary routines
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef AUX_HPP
#define AUX_HPP

#include <string>
#include <vector>
#include <iostream>
#include "matrix.hpp"

/** Determine reverse complement of a nucleic acid sequence */
std::string reverse_and_complement(const std::string &s);
/** Compute the sum of a matrix */
double norml0(const matrix_t &a);
/** Compute the L1-norm of a matrix */
double norml1(const matrix_t &a);
/** Compute the L2-norm of a matrix */
double norml2(const matrix_t &a);
/** Compute the Lp-norm of a matrix */
double normlp(const matrix_t &a, double p);
/** Compute the L1-norm of the difference of two matrices */
double norml1(const matrix_t &a, const matrix_t &b);
// /** Compute the L2-norm of the difference of two matrices */
// double norml2(const matrix_t &a, const matrix_t &b);
/** Compute the Lp-norm of the difference of two matrices */
double normlp(const matrix_t &a, const matrix_t &b, double p);

std::vector<size_t> shift_positions_left(
    const std::vector<size_t> original_inserts, size_t shift);
std::vector<size_t> shift_positions_right(
    const std::vector<size_t> original_inserts, size_t shift, size_t maxp);

/** Parse a comma separated list of ranges
 * Example: "1,2,5-7,2,5,30"  yields {1, 2, 5, 6, 7, 2, 5, 30} */
std::vector<size_t> parse_list(const std::string &s);

inline size_t ipow(size_t x, size_t y) {
  return static_cast<size_t>(pow(x, y));
}

/** A variant of getline that accepts both UNIX and windows line endings */
std::istream &safeGetline(std::istream &is, std::string &t);

/** Convert a range of character given by a pair of iterators to lower case */
template <typename Iter>
void range_tolower(Iter beg, Iter end) {
  for (Iter iter = beg; iter != end; ++iter) {
    *iter = std::tolower(*iter);
  }
}

/** Convert a string to lower case */
std::string string_tolower(const std::string &str);

/** Break long lines at whitespace to a maximal line length */
std::string limit_line_length(const std::string &x, size_t line_length);

/** Compute the logarithms of the sum of two logarithmic values */
inline double exp_add(double x, double y) {
  if (std::isinf(x) == -1 and std::isinf(y) == -1)
    return -std::numeric_limits<double>::infinity();
  double m = std::max(x, y);
  return log(exp(x - m) + exp(y - m)) + m;
}

/** Compute the logarithms of the difference of two logarithmic values */
inline double exp_diff(double x, double y) {
  if (std::isinf(x) == -1 and std::isinf(y) == -1)
    return -std::numeric_limits<double>::infinity();
  double m = std::max(x, y);
  return log(exp(x - m) - exp(y - m)) + m;
}

/** Parse a comma separated list of ranges
 * Example: "1,2,5-7,2,5,30"  yields {1, 2, 5, 6, 7, 2, 5, 30} */
std::vector<std::string> tokenize(const std::string &s,
                                  const std::string &delim);
std::string sha1hash(const std::string &s);

#endif
