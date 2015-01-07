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
 *       Filename:  matrix.hpp
 *
 *    Description:  Types and routines for Boost matrix and vectors classes
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using fp_t = double;
using index_t = int;

struct confusion_matrix {
  double true_positives;
  double false_negatives;
  double false_positives;
  double true_negatives;
};

confusion_matrix operator+(const confusion_matrix &m, double x);

using matrix_t = boost::numeric::ublas::matrix<fp_t>;
using vector_t = boost::numeric::ublas::vector<fp_t>;
using count_matrix_t = boost::numeric::ublas::matrix<size_t>;

using zero_matrix = boost::numeric::ublas::zero_matrix<fp_t>;
using identity_matrix = boost::numeric::ublas::identity_matrix<fp_t>;
using scalar_matrix = boost::numeric::ublas::scalar_matrix<fp_t>;

using zero_vector = boost::numeric::ublas::zero_vector<fp_t>;
using scalar_vector = boost::numeric::ublas::scalar_vector<fp_t>;

using zero_count_matrix = boost::numeric::ublas::zero_matrix<size_t>;
using identity_count_matrix = boost::numeric::ublas::identity_matrix<size_t>;
using scalar_count_matrix = boost::numeric::ublas::scalar_matrix<size_t>;

using count_vector_t = std::vector<size_t>;
std::string vec2string(const count_vector_t &v);

matrix_t operator+(const matrix_t &a, double x);
matrix_t operator-(const matrix_t &a, double x);

matrix_t read_emission(const std::string &path);

vector_t rowsums(const matrix_t &m);
vector_t colsums(const matrix_t &m);
matrix_t pum_consensus_pssm();
matrix_t qki_consensus_pssm();

void exp_transform(matrix_t &m);
void exp_transform(double &m);
void log_transform(matrix_t &m);
void log_transform(double &m);

void normalize(matrix_t &m);
void normalize(vector_t &m);

#endif
