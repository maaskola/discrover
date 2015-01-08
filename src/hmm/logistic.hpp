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
 *       Filename:  logistic.hpp
 *
 *    Description:  Code for various logistic calculations
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef LOGISTIC_HPP
#define LOGISTIC_HPP

/** The well-known sigmoid function, alternative calculation. */
inline double sigmoid_alternative(double x) {
  double y = exp(x);
  return y / (1 + y);
}

/** The well-known sigmoid function. */
inline double sigmoid(double x) { return 1 / (1 + exp(-x)); }

/** This gives 1 - sigmoid(x) */
inline double sigmoid_q(double x) { return 1 / (1 + exp(x)); }

inline double sigmoid_derivative_simple(double x) {
  double y = sigmoid(x);
  return y * (1 - y);
}

inline double sigmoid_derivative(double x) {
  return 1 / (exp(x) + 2 + exp(-x));
}

inline double logit(double p) { return log(p) - log(1 - p); }

inline double logistic_transformation(double x, double alpha, double beta) {
  return sigmoid(alpha * x + beta);
}

#endif
