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
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "aux.hpp"
#include "align.hpp"
#include "../timer.hpp"
#include "count.hpp"

using namespace std;

namespace Seeding {

  hash_map_t get_word_counts(const DataCollection &collection, size_t length, const options_t &options) {
    Timer t;

    size_t n_samples = 0;
    for(auto &series: collection)
      n_samples += series.sets.size();

    if(options.verbosity >= Verbosity::debug)
      cout << "Getting word counts for " << n_samples << " samples." << endl;

    hash_map_t counts;
    Stats::OccurrenceCounts default_stats(n_samples);
    for(auto &x: default_stats)
      x = 0;

    size_t idx = 0;
    for(auto &series: collection)
      for(auto &set: series)
        add_counts(set, length, counts, idx++, default_stats, options);

    double time = t.tock() * 1e-6;
    if(options.measure_runtime)
      cerr << "Getting word counts of length " << length << " took " << time << " seconds." << endl;

    return(counts);
  }

  size_t count_motif(const string &seq, const string &motif, const options_t &options) {
    size_t cnt = 0;
    auto qiter = begin(motif);
    auto qend = end(motif);
    auto riter = begin(seq);
    auto rend = end(seq);
    while((riter = search(riter, rend, qiter, qend, iupac_included)) != rend) {
      cnt++;
      if(not options.word_stats)
        return(1);
      riter++;
    }
    if(options.revcomp and (options.word_stats or cnt == 0)) {
      string rev_comp_motif = reverse_complement(motif);
      qiter = begin(rev_comp_motif);
      qend = end(rev_comp_motif);
      riter = begin(seq);
      while((riter = search(riter, rend, qiter, qend, iupac_included)) != rend) {
        cnt++;
        if(not options.word_stats)
          return(1);
        riter++;
      }
    }
    return(cnt);
  }

  /*
  Stats::OccurrenceCounts count_motif(const DataSeries &data_series, const string &motif, const options_t &options) {
    Stats::OccurrenceCounts stats(data_series.sets.size());
    for(auto &x: stats)
      x = 0;

    size_t idx = 0;
    for(auto &data_set: data_series) {
      for(auto &seq: data_set)
        stats(idx) += count_motif(seq.sequence, motif, options);
      idx++;
    }
    return(stats);
  }
  */

  Stats::OccurrenceCounts count_motif(const DataCollection &collection, const string &motif, const options_t &options) {
    size_t n_samples = 0;
    for(auto &series: collection)
      n_samples += series.sets.size();
    Stats::OccurrenceCounts stats(n_samples);
    for(auto &x: stats)
      x = 0;

    size_t idx = 0;
    for(auto &series: collection)
      for(auto &set: series) {
        for(auto &seq: set)
          stats(idx) += count_motif(seq.sequence, motif, options);
        idx++;
      }
    return(stats);
  }


  void print_counts(const hash_map_t &counts) {
    for(auto &iter: counts) {
      cout << iter.first;
      for(auto &x: iter.second)
        cout << "\t" << x;
      cout << endl;
    }
  }

  void add_counts(const string &seq_, size_t length, hash_map_t &counts, size_t idx, const Stats::OccurrenceCounts &default_stats, const options_t &options) {
    if(options.verbosity >= Verbosity::debug)
      cout << "Adding counts for sequence " << seq_ << endl;
    if(seq_.size() < length)
      return;

    vector<string> set;

    string seq = string_tolower(seq_);
    for(size_t i = 0; i + length <= seq.size(); i++) {
      string s = seq.substr(i, length);
      if(s.find_first_not_of("acgt") == string::npos) {
        if(options.revcomp) {
          string rc = reverse_complement(s);
          if(options.word_stats){
            set.push_back(s);
            set.push_back(rc);
          } else {
            if(lexicographical_compare(begin(s), end(s), begin(rc), end(rc)))
              set.push_back(rc);
            else
              set.push_back(s);
          }
        } else
          set.push_back(s);
      }
    }
    if(not options.word_stats) {
      sort(begin(set), end(set));
      set.resize(unique(begin(set), end(set)) - begin(set));
    }
    for(auto &w: set) {
      hash_map_t::iterator iter = counts.find(w);
      if(iter == end(counts)) {
        auto inserted = counts.insert({w,default_stats});
        iter = inserted.first;
      }
      iter->second(idx)++;
    }
  }

  void add_counts(const DataSet &data, size_t len, hash_map_t &counts, size_t idx, const Stats::OccurrenceCounts &default_stats, const options_t &options) {
    for(auto &seq: data)
      add_counts(seq.sequence, len, counts, idx, default_stats, options);
  }
}

