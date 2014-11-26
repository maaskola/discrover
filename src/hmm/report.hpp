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
 *       Filename:  report.hpp
 *
 *    Description:  Output routines
 *
 *        Created:  05/30/2012 06:42:25 PM
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef REPORT_HPP
#define REPORT_HPP

#include <iostream>
#include "hmm.hpp"
#include "../timer.hpp"

using namespace std;

namespace Evaluation {
  void print_table(ostream &ofs, const matrix_t m, const Data::Contrast &contrast, size_t width, size_t prec);

  void print_posterior(std::ostream &os, const HMM &hmm, const Data::Seq &seq);
  void count_report(ostream &ofs, const matrix_t counts, size_t motif_len, const Data::Contrast &contrast, double pseudo_count, bool limit_logp, const string &name, const string &prefix);

  void eval_contrast(const HMM &hmm, const Data::Contrast &contrast, ostream &ofs, bool limit_logp, const std::string &tag);

  /** Expected and Viterbi counts of occurrences and sites with the motifs as key*/
  struct ResultsCounts {
    typedef std::map<size_t, size_t> map_t;
    typedef std::map<size_t, double> float_map_t;
    float_map_t exp_sites;
    float_map_t exp_motifs;
    map_t viterbi_sites;
    map_t viterbi_motifs;
  };

  struct Result {
    struct Files {
      std::string summary;
      std::string viterbi;
      std::string table;
    };
    Files files;
  };

  /** Evaluate a single data set.
   * @return expected and Viterbi counts of occurrences and sites of all motifs
   */
  ResultsCounts evaluate_hmm_single_data_set(const HMM &hmm_,
      const Data::Set &dataset,
      ostream &out,
      ostream &v_out,
      ostream &occurrence_out,
      const Options::HMM &options);

  Result evaluate_hmm(const HMM &hmm_,
      const Data::Collection &collection,
      const std::string &tag,
      const Training::Tasks &tasks,
      const Options::HMM &options);
};
 
#endif

