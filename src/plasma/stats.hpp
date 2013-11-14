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
 *       Filename:  stats.hpp
 *
 *    Description:  Data structure to store counts
 *
 *        Created:  Thu Jun 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef  STATS_HPP
#define  STATS_HPP

#include "../matrix.hpp"

namespace Seeding {
  namespace Stats {
    typedef count_vector_t OccurrenceCounts;
    typedef matrix_t OccurrenceTable;
  }
}

#endif   /* ----- #ifndef STATS_HPP ----- */

