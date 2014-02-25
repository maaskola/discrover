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
 *       Filename:  association.hpp
 *
 *    Description:  Header for routines to calculate various measures of association
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef ASSOCIATION_HPP
#define ASSOCIATION_HPP

#include <vector>
#include "../matrix.hpp"

/** Computes mutual information of a contingency table
 * a,b,c,d are assumed to be non-negative
 */
double calc_mutual_information(double A, double B, double C, double D, bool normalize);
double calc_mutual_information(const confusion_matrix &m, bool normalize);
double calc_mutual_information(const matrix_t &matrix, double pseudo_count, bool normalize, bool correction, bool variance);

double calc_g_test_from_mi(double mi, double n);
double calc_g_test(const matrix_t &matrix, double pseudo_count=1.0);
double calc_log_likelihood_ratio(const matrix_t &matrix, double pseudo_count=1.0);

/** Computes mutual information of a contingency table
 * a,b,c,d are assumed to be non-negative
 */
double calc_matthews_correlation_coefficient(double tp, double fp, double fn, double tn);
double calc_matthews_correlation_coefficient(const confusion_matrix &m);

double calc_rank_information(vector_t posterior, double pseudo_count);

//TODO add: Cramer's V, Phi coefficient MCC, chi-square, DIPS

#endif

