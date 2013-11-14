/*
 * =====================================================================================
 *
 *       Filename:  mask.cpp
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

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "../aux.hpp"
#include "code.hpp"
#include "score.hpp"
#include "count.hpp"
#include "../timer.hpp"
#include "align.hpp"
// #include "GitSHA1.hpp"

using namespace std;

namespace Seeding {
  void remove_seqs_with_motif(const string &motif, DataSet &data_set, const Options &options) {
    vector<size_t> to_be_deleted;

    auto iter = begin(data_set.sequences);
    while(iter != end(data_set.sequences)) {
      if(search(begin(iter->sequence), end(iter->sequence), begin(motif), end(motif), iupac_included) != end(iter->sequence))
        to_be_deleted.push_back(distance(begin(data_set.sequences), iter));
      iter++;
    }
    if(options.verbosity >= Verbosity::verbose)
      cout << "Removing " << to_be_deleted.size() << " sequences that contain the motif." << endl;
    reverse(begin(to_be_deleted), end(to_be_deleted));
    for(auto &x: to_be_deleted) {
      data_set.set_size -= 1;
      data_set.seq_size -= data_set.sequences[x].sequence.size();
      if(options.verbosity >= Verbosity::debug)
        cout << "removing sequence " << x << " " << data_set.sequences[x].definition << " " << data_set.sequences[x].sequence << endl; //" " << data_set.sequences[x] << endl;
      data_set.sequences.erase(begin(data_set.sequences) + x);
    }
  }

  bool mask_motif_occurrences(const string &motif, string &seq, const Options &options, char mask_symbol) {
    if(options.verbosity >= Verbosity::debug)
      cout << "Check for masking of sequence " << seq << endl;
    set<size_t> occurrences;
    auto s_iter = begin(seq);
    while((s_iter = search(s_iter, end(seq), begin(motif), end(motif), iupac_included)) != end(seq))
      occurrences.insert(distance(begin(seq), s_iter++));
    if(options.revcomp) {
      string rc = reverse_complement(motif);
      s_iter = begin(seq);
      while((s_iter = search(s_iter, end(seq), begin(rc), end(rc), iupac_included)) != end(seq))
        occurrences.insert(distance(begin(seq), s_iter++));
    }
    if(not occurrences.empty()) {
      if(options.verbosity >= Verbosity::debug)
        cout << "Masking sequence " << seq << endl;
      for(auto &p: occurrences)
        for(size_t i = 0; i < motif.length(); i++)
          seq[p+i] = mask_symbol;
      if(options.verbosity >= Verbosity::debug)
        cout << "Masked sequence2 " << seq << endl;
      return(true);
    }
    return(false);
  }

  bool mask_motif_occurrences(const string &motif, Fasta::Entry &seq, const Options &options, char mask_symbol) {
    return(mask_motif_occurrences(motif, seq.sequence, options, mask_symbol));
  };

  void mask_motif_occurrences(const string &motif, DataSet &data_set, const Options &options) {
    for(auto &seq: data_set) {
      mask_motif_occurrences(motif, seq, options, 'n');
    }
  }

  void mask_motif_occurrences(const string &motif, DataSeries &data_series, const Options &options) {
    for(auto &data_set: data_series)
      mask_motif_occurrences(motif, data_set, options);
  }

  void remove_seqs_with_motif(const string &motif, DataSeries &data_series, const Options &options) {
    for(auto &data_set: data_series)
      remove_seqs_with_motif(motif, data_set, options);
  }

  void apply_mask(DataCollection &d, const string &motif, const Options &options) {
    switch(options.occurrence_filter) {
      case OccurrenceFilter::RemoveSequences:
        if(options.verbosity >= Verbosity::verbose)
          cout << "Removing sequences with " << motif << " occurrences." << endl;
        for(auto &series: d)
          remove_seqs_with_motif(motif, series, options);
        break;
      case OccurrenceFilter::MaskOccurrences:
        if(options.verbosity >= Verbosity::verbose)
          cout << "Masking " << motif << " occurrences." << endl;
        for(auto &series: d)
          mask_motif_occurrences(motif, series, options);
        break;
    }
  }


};

