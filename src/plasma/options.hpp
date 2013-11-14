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
 *       Filename:  options.hpp
 *
 *    Description:  Data structure to store program options of Plasma
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef OPTIONS_HPP
#define OPTIONS_HPP


#include <vector>
#include <string>
#include "measure.hpp"
#include "plasma_stats.hpp"
#include "specification.hpp"
#include "../verbosity.hpp"


namespace Seeding {
  enum class OccurrenceFilter {
    RemoveSequences,
    MaskOccurrences
  };
  std::istream &operator>>(std::istream &in, OccurrenceFilter &filter);
  std::ostream &operator<<(std::ostream &os, const OccurrenceFilter &filter);

  typedef Specification::Objective<Measures::Discrete::Measure> Objective;
  typedef std::vector<Objective> Objectives;

  std::ostream &operator<<(std::ostream &out, const Objectives &objectives);

  struct options_t {
    options_t();

    Specification::DataSets paths;
    Specification::Motifs motif_specifications;
    Objectives objectives;
    size_t n_threads;
    bool revcomp;
    bool strict;
    double pseudo_count;
    size_t n_seq;
    bool word_stats;
    bool measure_runtime;
    size_t n_motifs;
    OccurrenceFilter occurrence_filter;
    size_t max_candidates;
    std::vector<size_t> degeneracies;
    double rel_degeneracy;
    bool per_degeneracy;
    bool keep_all;
    Verbosity verbosity;
    bool dump_viterbi;
    bool no_enrichment_filter;
  };
}


#endif   /* ----- #ifndef OPTIONS_HPP  ----- */

