/* =====================================================================================
 * Copyright (c) 2015, Jonas Maaskola
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
 *       Filename:  random_distributions.cpp
 *
 *    Description:  Shared random distriubtions
 *
 *        Created:  Thu Jan 15 06:44:43 2015 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef RANDOM_DISTRIBUTIONS_HPP
#define RANDOM_DISTRIBUTIONS_HPP

#include <random>

struct RandomDistribution {
  static std::uniform_int_distribution<size_t> Uniform;
  static std::uniform_int_distribution<size_t> Binary;
  static std::uniform_int_distribution<size_t> Nucleotide;
  static std::uniform_real_distribution<double> Probability;
};

#endif
