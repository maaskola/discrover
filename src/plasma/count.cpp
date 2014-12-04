/*
 * =====================================================================================
 *
 *       Filename:  count.cpp
 *
 *    Description:  
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "../aux.hpp"
#include "code.hpp"
#include "data.hpp"
#include "../timer.hpp"
#include "count.hpp"

using namespace std;

namespace Seeding {

  hash_map_t get_word_counts(const Collection &collection, size_t length, const Options &options) {
    Timer t;

    size_t n_samples = 0;
    for(auto &contrast: collection)
      n_samples += contrast.sets.size();

    if(options.verbosity >= Verbosity::debug)
      cout << "Getting word counts for " << n_samples << " samples." << endl;

    hash_map_t counts;
    count_vector_t default_stats(n_samples);
    for(auto &x: default_stats)
      x = 0;

    size_t idx = 0;
    for(auto &contrast: collection)
      for(auto &dataset: contrast)
        add_counts(dataset, length, counts, idx++, default_stats, options);

    double time = t.tock() * 1e-6;
    if(options.measure_runtime)
      cerr << "Getting word counts of length " << length << " took " << time << " seconds." << endl;

    return(counts);
  }

  size_t count_motif(const string &seq, const string &motif, const Options &options) {
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
  count_vector_t count_motif(const Contrast &contrast, const string &motif, const Options &options) {
    count_vector_t stats(contrast.sets.size());
    for(auto &x: stats)
      x = 0;

    size_t idx = 0;
    for(auto &dataset: contrast) {
      for(auto &seq: dataset)
        stats(idx) += count_motif(seq.sequence, motif, options);
      idx++;
    }
    return(stats);
  }
  */

  count_vector_t count_motif(const Collection &collection, const string &motif, const Options &options) {
    size_t n_samples = 0;
    for(auto &contrast: collection)
      n_samples += contrast.sets.size();
    count_vector_t stats(n_samples);
    for(auto &x: stats)
      x = 0;

    size_t idx = 0;
    for(auto &contrast: collection)
      for(auto &dataset: contrast) {
        for(auto &seq: dataset)
          stats[idx] += count_motif(seq.sequence, motif, options);
        idx++;
      }
    return(stats);
  }


  void print_counts(const hash_map_t &counts) {
    for(auto &iter: counts) {
      cout << decode(iter.first);
      for(auto &x: iter.second)
        cout << "\t" << x;
      cout << endl;
    }
  }

  void add_counts(const string &seq_, size_t length, hash_map_t &counts, size_t idx, const count_vector_t &default_stats, const Options &options) {
    if(options.verbosity >= Verbosity::debug)
      cout << "Adding counts for sequence " << seq_ << endl;
    if(seq_.size() < length)
      return;

    vector<hash_map_t::key_type> set;
    vector<hash_map_t::key_type> set2;

    Timer t1, t2;
    double time1, time2;
    if (true) {
      t1.tick();
      seq_type seq = encode(seq_);

      auto b1 = begin(seq);
      auto b2 = b1;
      std::advance(b2, length);
      const auto e = end(seq);

      do {
        seq_type s(b1, b2);
        if (find_if(b1, b2, degenerate_nucleotide) == b2) {
          if (options.revcomp) {
            auto rc = iupac_reverse_complement(s);
            if (options.word_stats) {
              set2.push_back(s);
              set2.push_back(rc);
            } else {
              if (lexicographical_compare(begin(s), end(s), begin(rc), end(rc)))
                set2.push_back(rc);
              else
                set2.push_back(s);
            }
          } else
            set2.push_back(s);
        }
        ++b1;
      } while (b2++ != e);
      time1 = t1.tock();
    }

    t2.tick();
    string seq = string_tolower(seq_);
    for(size_t i = 0; i + length <= seq.size(); i++) {
      string s = seq.substr(i, length);
      if(s.find_first_not_of("acgt") == string::npos) {
        if(options.revcomp) {
          string rc = reverse_complement(s);
          if(options.word_stats){
            set.push_back(encode(s));
            set.push_back(encode(rc));
          } else {
            if(lexicographical_compare(begin(s), end(s), begin(rc), end(rc)))
              set.push_back(encode(rc));
            else
              set.push_back(encode(s));
          }
        } else
          set.push_back(encode(s));
      }
    }
    time2 = t2.tock();

    if(options.verbosity >= Verbosity::debug) {
      cerr << "time1 = " << time1 << " time2 = " << time2 << std::endl;
      cout << "set.size() = " << set.size() << endl;
      cout << "set2.size() = " << set2.size() << endl;
      auto x = begin(set);
      auto y = begin(set2);
      for( ; x != end(set) and y != end(set2); x++, y++)
        cout << decode(*x) << "\t" << decode(*y) << endl;
      for( ; x != end(set); x++)
        cout << decode(*x) << "\t" << "-X-" << endl;
      for( ; y != end(set2); y++)
        cout << "-X-" << "\t" << decode(*y) << endl;
    }
    assert(set.size() == set2.size());
    for(auto x = begin(set), y = begin(set2); x != end(set) and y != end(set2); ++x, ++y)
      assert(*x == *y);

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
      iter->second[idx]++;
    }
  }

  void add_counts(const Set &dataset, size_t len, hash_map_t &counts, size_t idx, const count_vector_t &default_stats, const Options &options) {
    for(auto &seq: dataset)
      add_counts(seq.sequence, len, counts, idx, default_stats, options);
  }
}
