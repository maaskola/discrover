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

#include "random_distributions.hpp"

using namespace std;

uniform_int_distribution<size_t> RandomDistribution::Uniform;
uniform_int_distribution<size_t> RandomDistribution::Binary(0, 1);
uniform_int_distribution<size_t> RandomDistribution::Nucleotide(0, 3);
uniform_real_distribution<double> RandomDistribution::Probability(0, 1);
