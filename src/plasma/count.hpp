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
 *       Filename:  count.hpp
 *
 *    Description:  Word counting routines
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef COUNT_HPP
#define COUNT_HPP

#include <iostream>
#include <algorithm>
#include "options.hpp"

namespace Seeding {

  /** Count occurrences of non-degenerate words; i.e. of words of the given length, consisting only of ACGT */
  void add_counts(const std::string &seq, size_t length, hash_map_t &counts, size_t idx, const count_vector_t &default_stats, const Options &options);
  void add_counts(const DataSet &data, size_t len, hash_map_t &counts, size_t idx, const count_vector_t &default_stats, const Options &options);

  typedef uint16_t seq_idx_t;
  typedef std::map<std::string, std::list<seq_idx_t>> Index;

  hash_map_t get_word_counts(const DataCollection &data_series, size_t length, const Options &options);
  // count_vector_t count_motif(const DataSeries &data_series, const std::string &motif, const Options &options);
  count_vector_t count_motif(const DataCollection &collection, const std::string &motif, const Options &options);

  void print_counts(const hash_map_t &counts);

}

#endif   /* ----- #ifndef COUNT_HPP  ----- */

