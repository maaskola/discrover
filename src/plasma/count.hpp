/*
 * =====================================================================================
 *
 *       Filename:  count.cpp
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


#ifndef  COUNT_HPP
#define  COUNT_HPP

#include <iostream>
#include <algorithm>
#include "options.hpp"

namespace Plasma {

  /** Count occurrences of non-degenerate words; i.e. of words of the given length, consisting only of ACGT */
  void add_counts(const std::string &seq, size_t length, hash_map_t &counts, size_t idx, const Stats::OccurrenceCounts &default_stats, const options_t &options);
  void add_counts(const DataSet &data, size_t len, hash_map_t &counts, size_t idx, const Stats::OccurrenceCounts &default_stats, const options_t &options);

  typedef uint16_t seq_idx_t;
  typedef std::map<std::string, std::list<seq_idx_t>> Index;

  hash_map_t get_word_counts(const DataCollection &data_series, size_t length, const options_t &options);
  // Stats::OccurrenceCounts count_motif(const DataSeries &data_series, const std::string &motif, const options_t &options);
  Stats::OccurrenceCounts count_motif(const DataCollection &collection, const std::string &motif, const options_t &options);

  void print_counts(const hash_map_t &counts);

}

#endif   /* ----- #ifndef COUNT_HPP  ----- */

