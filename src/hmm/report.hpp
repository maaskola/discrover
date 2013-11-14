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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef REPORT_HPP
#define REPORT_HPP

#include <iostream>
#include "hmm.hpp"
#include "../timer.hpp"

using namespace std;

void print_table(ostream &ofs, const matrix_t m, const Data::Series &data, size_t width, size_t prec);

void count_report(ostream &ofs, const matrix_t counts, size_t motif_len, const Data::Series &data, double pseudo_count, bool limit_logp, const string &name, const string &prefix);

void eval_contrast(const HMM &hmm, const Data::Series &data, ostream &ofs, bool limit_logp, const std::string &tag);

/** Expected and Viterbi counts of occurrences and sites with the motifs as key*/
struct ResultsCounts {
  typedef std::map<size_t, size_t> map_t;
  typedef std::map<size_t, double> float_map_t;
  float_map_t exp_sites;
  float_map_t exp_motifs;
  map_t viterbi_sites;
  map_t viterbi_motifs;
};

/** Evaluate a single data set.
 * @return expected and Viterbi counts of occurrences and sites of all motifs
 */
ResultsCounts evaluate_hmm_single_data_set(HMM &hmm_,
    const Data::Set &data,
    ostream &out,
    ostream &v_out,
    ostream &occurrence_out,
    const hmm_options &options);

void evaluate_hmm(HMM &hmm_,
    const Data::Collection &all_data,
    const Data::Collection &training_data,
    const Data::Collection &test_data,
    const Training::Tasks &tasks,
    const hmm_options &options);
 
#endif

