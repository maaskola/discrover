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
 *       Filename:  plasma_stats.hpp
 *
 *    Description:  Data structures to store count data
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef  PLASMA_STATS_HPP
#define  PLASMA_STATS_HPP

#include <unordered_map>
#include <map>
#include <string>
#include "align.hpp"
#include "../matrix.hpp"

namespace std {

template <>
struct hash<seq_type> {
  std::size_t operator()(const seq_type& v) const {
    using std::size_t;
    using std::hash;
    using std::string;

    // slow but compatible with previous, string-based implementation
    return hash<string>()(decode(v));

    // faster
    // size_t h = 0;
    // for(auto &x: v)
    //   h = h ^ x;
    // return h;

    // faster, likely the same as above, if the implementation uses identity hashes for integral types
    // size_t h = 0;
    // for(auto &x: v)
    //   h = h ^ hash<seq_type::value_type>()(x);
    // return h;

  }
};
}

namespace Seeding {
  using hash_map_t = std::unordered_map<seq_type, count_vector_t>;
  using score_map_t = std::unordered_map<seq_type, double>;
  using rev_map_t = std::multimap<double, seq_type, std::greater<double>>;
}

#endif   /* ----- #ifndef PLASMA_STATS_HPP ----- */
