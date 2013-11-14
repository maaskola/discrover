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
 *       Filename:  find.hpp
 *
 *    Description:  Routines to find IUPAC regular expression type motifs
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef FIND_HPP
#define FIND_HPP

#include "options.hpp"
#include "results.hpp"
#include "align.hpp"

namespace Seeding {
  struct Plasma {
    options_t options;
    DataCollection collection;
    bool index_ready;
    bool needs_rebuilding;
    NucleotideIndex<size_t,size_t> index;
    Plasma(const options_t &options_t);
    Plasma(const DataCollection &collection_, const options_t &opt);
    Results find_breadth(size_t length, const Objective &objective);
    Results find_all(const Specification::Motif &motif, const Objective &objective, size_t n_motifs) const;
    Results find_multiple(const Specification::Motif &motif, const Objective &objective, size_t n_motifs) const;
    void apply_mask(const std::string &motif);
    void apply_mask(const Result &result);
    void apply_mask(const Results &results);
    void rebuild_index();
    Results find(const Specification::Motif &motif, const Objectives &objectives, bool doreport=true) const;
  };

  void report(std::ostream &os, const Objective &objective, const std::string &motif, const DataCollection &collection, const options_t &options);
  void report(std::ostream &os, const Result &result, const DataCollection &collection, const options_t &options);
  void viterbi_dump(const std::string &motif, const DataCollection &collection, std::ostream &out, const options_t &options);
}

#endif   /* ----- #ifndef FIND_HPP  ----- */

