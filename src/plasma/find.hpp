/*
 * =====================================================================================
 *
 *       Filename:  find.hpp
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


#ifndef FIND_HPP
#define FIND_HPP

#include "options.hpp"
#include "results.hpp"
#include "align.hpp"

namespace Plasma {
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

