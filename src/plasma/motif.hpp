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
 *       Filename:  motif.hpp
 *
 *    Description:  Data structure for regular expression type motifs
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef MOTIF_HPP
#define MOTIF_HPP

#include <cstdint>
#include <vector>
#include <string>
#include "code.hpp"
#include "../verbosity.hpp"

namespace Seeding {
  double information_content(const std::string &motif);
  size_t motif_degeneracy(const std::string &motif);

  std::vector<seq_type> all_generalizations(const seq_type &motif);
};

#endif   /* ----- #ifndef MOTIF_HPP ----- */

