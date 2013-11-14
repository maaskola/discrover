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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/algorithm/string.hpp>
#include "suffix.hpp"

#include "../aux.hpp"
#include "data.hpp"

bool iupac_included(char q, char r);

std::string iupac_reverse_complement(const std::string &s);

std::string read_fasta_with_boundaries(const std::vector<std::string> &paths, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set, size_t n_seq=0);

template <class data_t, class idx_t=size_t, class lcp_t=size_t>
void list_occurrences(const data_t &x, const data_t &q) {
  // generate suffix array, LCP, and JMP tables
  std::vector<idx_t> sa = gen_suffix_array<idx_t>(begin(x), end(x));
  std::vector<lcp_t> lcp = gen_lcp<lcp_t>(begin(x), end(x), sa);
  std::vector<idx_t> jmp = gen_jmp<idx_t>(lcp);

  std::list<idx_t> hits = match(begin(q), end(q), begin(x), end(x), sa, lcp, jmp, iupac_included);
  for(auto &hit: hits) {
    std::cout << "Found match at position " << hit << ":\t";
    for(size_t i = 0; i < q.size(); i++)
      std::cout << x[hit+i];
    std::cout << std::endl;

  }
  if(hits.empty())
    std::cout << "No hits found." << std::endl;
}

template <class data_t, class idx_t=size_t, class lcp_t=size_t>
class Index {
  public:
    Index(const data_t &x, Verbosity verbosity) :
      data(x),
      sa(gen_suffix_array<idx_t>(begin(data), end(data), verbosity)),
      lcp(gen_lcp<lcp_t>(begin(data), end(data), sa, verbosity)),
      jmp(gen_jmp<idx_t>(lcp, verbosity)) {
        if(verbosity >= Verbosity::verbose) {
          std::cout << "SA hash = " << boost::hash_range(begin(sa), end(sa)) << std::endl;
          std::cout << "LCP hash = " << boost::hash_range(begin(lcp), end(lcp)) << std::endl;
          std::cout << "JMP hash = " << boost::hash_range(begin(jmp), end(jmp)) << std::endl;
        }
      };
    template<class Cmp=std::equal_to<typename data_t::value_type>> std::list<idx_t> find_matches(const data_t &query, Cmp cmp=Cmp()) const {
      return(match(begin(query), end(query), begin(data), end(data), sa, lcp, jmp, cmp));
    };
  private:
    data_t data;               // the original data
    std::vector<idx_t> sa;     // suffix array
    std::vector<lcp_t> lcp;    // LCP table
    std::vector<idx_t> jmp;    // JMP table
};

template <class data_t, class idx_t=size_t, class lcp_t=size_t>
class BidirectionalIndex{
  public:
    BidirectionalIndex(const data_t &x) :
      forward_index(x),
      backward_index(data_t(x.rbegin(), x.rend())) {
      };
    template<class Cmp=std::equal_to<typename data_t::value_type>> std::list<idx_t> find_matches(const data_t &query, Cmp cmp=Cmp()) const {
      return(forward_index.find_matches(query, cmp));
    };
  private:
    Index<data_t, idx_t, lcp_t> forward_index, backward_index;
};

// std::string collapse_data_series(const Seeding::DataSeries &data_series, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set);
std::string collapse_data_collection(const Seeding::DataCollection &collection, std::vector<size_t> &pos2seq, std::vector<size_t> &seq2set, std::vector<size_t> &set2series);

template <class idx_t=size_t, class lcp_t=size_t, class index_t=Index<std::string, idx_t, lcp_t>>
class NucleotideIndex {
  struct Hit {
    size_t file_idx;
    size_t entry_idx;
    size_t pos;
  };
  public:
    NucleotideIndex(const NucleotideIndex &i) :
      paths(i.paths),
      pos2seq(i.pos2seq),
      seq2set(i.seq2set),
      set2series(i.set2series),
      index(i.index) { };
    NucleotideIndex(Verbosity verbosity=Verbosity::info) :
      paths(),
      pos2seq(),
      seq2set(),
      set2series(),
      index("", verbosity) { };
    NucleotideIndex(const Seeding::DataCollection &collection, Verbosity verbosity) :
      paths(),
      pos2seq(),
      seq2set(),
      set2series(),
      index(collapse_data_collection(collection, pos2seq, seq2set, set2series), verbosity) {
        for(auto &series: collection)
          for(auto &set: series)
            paths.push_back(set.path);
      };
//    NucleotideIndex(const Seeding::DataSeries &data_series, Verbosity verbosity) :
//      paths(),
//      pos2seq(),
//      seq2set(),
//      index(collapse_data_series(data_series, pos2seq, seq2set), verbosity) {
//        for(auto &data_set: data_series)
//          paths.push_back(data_set.path);
//      };
    NucleotideIndex(const std::vector<std::string> &paths_, size_t nseq=0) :
      paths(paths_),
      pos2seq(),
      seq2set(),
      set2series(),
      index(read_fasta_with_boundaries(paths, pos2seq, seq2set, nseq)) {
      };
    std::vector<size_t> word_hits_by_file(const std::string &query) const {
      std::vector<size_t> counts(paths.size(),0);
      for(auto &p: index.find_matches(query, iupac_included)) {
        size_t seqIdx = pos2seq[p];
        size_t fileIdx = seq2set[seqIdx];
        counts[fileIdx]++;
      }
      return(counts);
    };
    std::vector<size_t> seq_hits_by_file(const std::string &query, bool revcomp=false) const {
      std::vector<size_t> seqs;
      for(auto &p: index.find_matches(query, iupac_included))
        seqs.push_back(pos2seq[p]);
      if(revcomp)
        for(auto &p: index.find_matches(iupac_reverse_complement(query), iupac_included))
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
    std::vector<size_t> pos2seq, seq2set, set2series;
    index_t index;
};
template <class idx_t=size_t, class lcp_t=size_t>
using BidirectionalNucleotideIndex = NucleotideIndex<idx_t, lcp_t, BidirectionalIndex<std::string, idx_t, lcp_t>>;

#endif

