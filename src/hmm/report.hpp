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

class Evaluator {
  HMM hmm;

public:
  Evaluator(const HMM &hmm);

  struct Result {
    struct Files {
      std::string summary;
      std::string viterbi;
      std::string bed;
      std::string table;
      std::vector<std::string> logos;
    };
    Files files;
  };

  Result report(const Data::Collection &collection, const std::string &tag,
                const Training::Tasks &tasks,
                const Options::HMM &options) const;

private:
  void print_posterior(std::ostream &os, const vector_t &scale,
                       const matrix_t &f, const matrix_t &b) const;
  void eval_contrast(std::ostream &ofs, const Data::Contrast &contrast,
                     bool limit_logp, const std::string &tag) const;

  /** Expected and Viterbi counts of occurrences and sites with motifs as key
   **/
  struct ResultsCounts {
    using map_t = std::map<size_t, size_t>;
    using float_map_t = std::map<size_t, double>;
    float_map_t exp_sites;
    float_map_t exp_motifs;
    map_t viterbi_sites;
    map_t viterbi_motifs;
  };

  /** Evaluate a single data set.
   * @return expected and Viterbi counts of occurrences and sites of all motifs
   */
  ResultsCounts evaluate_dataset(const Data::Set &dataset, std::ostream &out,
                                 std::ostream &v_out, std::ostream &occ_out,
                                 std::ostream &motif_out, std::ostream &bed_out,
                                 const Options::HMM &options) const;

#if CAIRO_FOUND
  std::vector<std::string> generate_logos(const std::string &path_stem,
                                          const Options::HMM &options) const;
#endif
};

#endif
