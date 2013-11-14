/* =====================================================================================
 * Copyright (c) 2012, Jonas Maaskola
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
 *       Filename:  mic_impl.hpp
 *
 *    Description:  Routines to compute the MIC
 *
 *        Created:  05/30/2012 06:42:25 PM
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef MIC_IMPL_HPP
#define MIC_IMPL_HPP

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include "association.hpp"
#include "../verbosity.hpp"

namespace MIC {

  typedef double score_t;
  //typedef size_t data_t;
  typedef double data_t;
  typedef std::vector<data_t> sample_t;
  typedef std::vector<size_t> bounds_t;
  typedef sample_t stats_t;

  bounds_t sample(const sample_t &d, const std::set<size_t> &fixed_bounds, size_t n, size_t n_iter, Verbosity verbosity);

};

#endif

