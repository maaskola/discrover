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
 *       Filename:  statistics.hpp
 *
 *    Description:  Some code for random number generation using Boost routines
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef STATISTICS_HPP
#define STATISTICS_HPP


#include <boost/random.hpp>

double rgamma( double mean, double variance, boost::mt19937& rng );
std::vector<double> rdirichlet(std::vector<double> alpha, boost::mt19937& rng);
std::vector<double> runiform(size_t n, boost::mt19937& rng);

//  boost::mt19937 rng;
//  for(size_t i = 0; i < n; i++)
//    cout << runiform(k, rng) << endl;

#endif     /* -----  STATISTICS_HPP  ----- */

