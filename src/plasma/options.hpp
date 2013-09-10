/*
 * =====================================================================================
 *
 *       Filename:  options.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
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


namespace Plasma {
  enum class OccurrenceFilter {
    remove_seq,
    mask_occurrence
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

