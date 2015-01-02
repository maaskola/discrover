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
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <vector>
#include <string>
#include "measure.hpp"
#include "specification.hpp"
#include "../verbosity.hpp"

namespace Seeding {
  enum class OccurrenceFilter {
    RemoveSequences,
    MaskOccurrences
  };

  std::istream &operator>>(std::istream &in, OccurrenceFilter &filter);
  std::ostream &operator<<(std::ostream &os, const OccurrenceFilter &filter);

  enum class Algorithm {
    Plasma = (1u << 1),
    ExternalDREME = (1u << 2),
    MCMC = (1u << 3)
  };

  inline Algorithm operator|(Algorithm a, Algorithm b)
  {return static_cast<Algorithm>(static_cast<int>(a) | static_cast<int>(b));}
  inline Algorithm operator&(Algorithm a, Algorithm b)
  {return static_cast<Algorithm>(static_cast<int>(a) & static_cast<int>(b));}

  std::istream &operator>>(std::istream &in, Algorithm &algorithm);
  std::ostream &operator<<(std::ostream &os, const Algorithm &algorithm);

  using Objective = Specification::Objective<Measures::Discrete::Measure>;
  using Objectives = std::vector<Objective>;

  Objective objective_for_motif(const Objectives &objectives, const Specification::Motif &motif);

  std::ostream &operator<<(std::ostream &out, const Objective &objective);

  struct Options {
    struct Plasma {
      Plasma();
      size_t max_candidates;
      std::vector<size_t> degeneracies;
      double rel_degeneracy;
      bool per_degeneracy;
    };
    struct MCMC {
      MCMC();
      size_t max_iter;
      double temperature;
      size_t n_parallel;
      unsigned int random_salt;
    };


    Options();

    Specification::Sets paths;
    Specification::Motifs motif_specifications;
    Objectives objectives;

    Algorithm algorithm;

    Plasma plasma;
    MCMC mcmc;

    size_t n_threads;
    bool revcomp;
    bool strict;
    double pseudo_count;
    bool weighting;
    size_t n_seq;
    bool word_stats;
    bool measure_runtime;
    OccurrenceFilter occurrence_filter;
    bool only_best;
    Verbosity verbosity;
    bool dump_viterbi;
    bool no_enrichment_filter;
    bool fixed_motif_space_mode;
    bool allow_iupac_wildcards;

    std::string label;
    bool pdf_logo;
    bool png_logo;
  };
}


#endif   /* ----- #ifndef OPTIONS_HPP  ----- */

