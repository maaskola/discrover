/*
 * =====================================================================================
 *
 *       Filename:  mask.hpp
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


#ifndef  MASK_HPP
#define  MASK_HPP

#include "options.hpp"

namespace Plasma {
  void remove_seqs_with_motif(const std::string &motif, DataSet &data_set, const options_t &options);
  bool mask_motif_occurrences(const std::string &motif, std::string &seq, const options_t &options, char mask_symbol='n');
  void mask_motif_occurrences(const std::string &motif, DataSet &data_set, const options_t &options);
  void mask_motif_occurrences(const std::string &motif, DataSeries &data_series, const options_t &options);
  void remove_seqs_with_motif(const std::string &motif, DataSeries &data_series, const options_t &options);
  void apply_mask(DataCollection &d, const std::string &motif, const options_t &options);
};

#endif   /* ----- #ifndef MASK_HPP  ----- */

