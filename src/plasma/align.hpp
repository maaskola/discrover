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
 *       Filename:  align.hpp
 *
 *    Description:  Alignment routines based on suffix arrays
 *
 *        Created:  Thu Aug 8 22:23:44 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <iostream>
#include <cstring>
#include <vector>
#include "suffix.hpp"
#include "data.hpp"

using symbol_t = uint8_t;
void add_sequence(std::vector<symbol_t> &s, const std::string &seq);

template <typename T> bool binary_and_not_null(T a, T b) {
  return (a & b) != 0;
}

std::string iupac_reverse_complement(const std::string &s);

std::string read_fasta_with_boundaries(const std::vector<std::string> &paths, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set, size_t n_seq=0);

template <class data_t, class idx_t=size_t, class lcp_t=size_t>
class Index {
  public:
    Index(const data_t &x, Verbosity verbosity) :
      data(x),
      // generate suffix array, LCP, and JMP tables
      sa(gen_suffix_array<idx_t>(begin(data), end(data), verbosity)),
      lcp(gen_lcp<lcp_t>(begin(data), end(data), sa, verbosity)),
      jmp(gen_jmp<idx_t>(lcp, verbosity)) { };

    template<class Cmp=std::equal_to<typename data_t::value_type>>
    std::vector<idx_t> find_matches(const data_t &query, Cmp cmp=Cmp()) const {
      return(match(begin(query), end(query), begin(data), end(data), sa, lcp, jmp, cmp));
    };
  private:
    data_t data;               // the original data
    std::vector<idx_t> sa;     // suffix array
    std::vector<lcp_t> lcp;    // LCP table
    std::vector<idx_t> jmp;    // JMP table
};

// std::string collapse_contrast(const Seeding::Contrast &contrast, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set);
std::vector<symbol_t> collapse_collection(const Seeding::Collection &collection, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set, std::vector<size_t> &set2contrast);

template <class idx_t=size_t, class lcp_t=size_t, class index_t=Index<std::vector<symbol_t>, idx_t, lcp_t>>
class NucleotideIndex {
  public:
    NucleotideIndex(const NucleotideIndex &i) :
      paths(i.paths),
      pos2seq(i.pos2seq),
      seq2set(i.seq2set),
      set2contrast(i.set2contrast),
      index(i.index) { };

    NucleotideIndex(Verbosity verbosity=Verbosity::info) :
      paths(),
      pos2seq(),
      seq2set(),
      set2contrast(),
      index({}, verbosity) { };

    NucleotideIndex(const Seeding::Collection &collection, Verbosity verbosity) :
      paths(),
      pos2seq(),
      seq2set(),
      set2contrast(),
      index(collapse_collection(collection, pos2seq, seq2set, set2contrast), verbosity) {
        for(auto &contrast: collection)
          for(auto &dataset: contrast)
            paths.push_back(dataset.path);
      };
    std::vector<size_t> word_hits_by_file(const std::string &query) const {
      std::vector<symbol_t> q;
      add_sequence(q, query);
      std::vector<size_t> counts(paths.size(),0);
      for(auto &p: index.find_matches(q, binary_and_not_null<symbol_t>)) {
        size_t seqIdx = pos2seq[p];
        size_t fileIdx = seq2set[seqIdx];
        counts[fileIdx]++;
      }
      return(counts);
    };
    std::vector<size_t> seq_hits_by_file(const std::string &query, bool revcomp=false) const {
      std::vector<symbol_t> q, qrc;
      add_sequence(q, query);
      add_sequence(qrc, iupac_reverse_complement(query));
      std::vector<size_t> seqs;
      for(auto &p: index.find_matches(q, binary_and_not_null<symbol_t>))
        seqs.push_back(pos2seq[p]);
      if(revcomp)
        for(auto &p: index.find_matches(qrc, binary_and_not_null<symbol_t>))
          seqs.push_back(pos2seq[p]);

      std::sort(begin(seqs), end(seqs));
      seqs.resize(std::unique(begin(seqs), end(seqs)) - begin(seqs));

      std::vector<size_t> counts(paths.size(),0);
      for(auto &s: seqs)
        counts[seq2set[s]]++;
      return(counts);
    };

  private:
    std::vector<std::string> paths;
    std::vector<size_t> pos2seq, seq2set, set2contrast;
    index_t index;
};

#endif
