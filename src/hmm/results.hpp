/* =====================================================================================
 * Copyright (c) 2013, Jonas Maaskola
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
 *       Filename:  results.hpp
 *
 *    Description:  Data structures to store Discrover results
 *
 *        Created:  Thu Mar 07 20:30:35 2014 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef HMM_RESULTS_HPP
#define HMM_RESULTS_HPP

#include <cstddef>
#include <vector>
#include <string>

namespace Training {

  struct State {
    State(size_t n=0);
    int center;
    std::vector<std::vector<double>> scores;
  };

  struct Result {
    Result();
    State state;
    double delta;
    std::string parameter_file;
    // bool success; // TODO add a flag to indicate training success
  };
}

#endif   /* ----- #ifndef RESULTS_HPP  ----- */

