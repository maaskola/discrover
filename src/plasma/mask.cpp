/*
 * =====================================================================================
 *
 *       Filename:  mask.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include "../aux.hpp"
#include "code.hpp"
#include "score.hpp"
#include "count.hpp"
#include "../timer.hpp"
#include "align.hpp"

using namespace std;

const char MASK_SYMBOL = '-';

const bool update_sizes_on_removal = true;

namespace Seeding {
void remove_seqs_with_motif(const string &motif, Set &dataset,
                            const Options &options) {
  auto pred = [&motif](const Fasta::Entry &entry) {
    return search(begin(entry.sequence), end(entry.sequence), begin(motif),
                  end(motif), iupac_included) != end(entry.sequence);
  };
  dataset.sequences.erase(
      remove_if(begin(dataset.sequences), end(dataset.sequences), pred),
      end(dataset.sequences));
  if (update_sizes_on_removal) {
    dataset.set_size = dataset.sequences.size();
    dataset.seq_size = 0;
    for (auto &entry : dataset.sequences)
      dataset.seq_size += entry.sequence.size();
  }
}

bool mask_motif_occurrences(const string &motif, string &seq,
                            const Options &options, char mask_symbol) {
  if (options.verbosity >= Verbosity::debug)
    cout << "Check for masking of sequence " << seq << endl;
  set<size_t> occurrences;
  auto s_iter = begin(seq);
  while ((s_iter = search(s_iter, end(seq), begin(motif), end(motif),
                          iupac_included)) != end(seq))
    occurrences.insert(distance(begin(seq), s_iter++));
  if (options.revcomp) {
    string rc = reverse_complement(motif);
    s_iter = begin(seq);
    while ((s_iter = search(s_iter, end(seq), begin(rc), end(rc),
                            iupac_included)) != end(seq))
      occurrences.insert(distance(begin(seq), s_iter++));
  }
  if (not occurrences.empty()) {
    if (options.verbosity >= Verbosity::debug)
      cout << "Masking sequence " << seq << endl;
    for (auto &p : occurrences)
      for (size_t i = 0; i < motif.length(); i++)
        seq[p + i] = mask_symbol;
    if (options.verbosity >= Verbosity::debug)
      cout << "Masked sequence2 " << seq << endl;
    return true;
  }
  return false;
}

bool mask_motif_occurrences(const string &motif, Fasta::Entry &seq,
                            const Options &options, char mask_symbol) {
  return mask_motif_occurrences(motif, seq.sequence, options, mask_symbol);
};

void mask_motif_occurrences(const string &motif, Set &dataset,
                            const Options &options) {
  for (auto &seq : dataset) {
    mask_motif_occurrences(motif, seq, options, MASK_SYMBOL);
  }
}

void mask_motif_occurrences(const string &motif, Contrast &contrast,
                            const Options &options) {
  for (auto &dataset : contrast)
    mask_motif_occurrences(motif, dataset, options);
}

void remove_seqs_with_motif(const string &motif, Contrast &contrast,
                            const Options &options) {
  for (auto &dataset : contrast)
    remove_seqs_with_motif(motif, dataset, options);
  if (update_sizes_on_removal) {
    contrast.seq_size = 0;
    contrast.set_size = 0;
    for (auto &dataset : contrast) {
      contrast.seq_size += dataset.seq_size;
      contrast.set_size += dataset.set_size;
    }
  }
}

void apply_mask(Collection &d, const string &motif, const Options &options) {
  switch (options.occurrence_filter) {
    case OccurrenceFilter::RemoveSequences:
      if (options.verbosity >= Verbosity::verbose)
        cout << "Removing sequences with " << motif << " occurrences." << endl;
      for (auto &contrast : d)
        remove_seqs_with_motif(motif, contrast, options);
      if (update_sizes_on_removal) {
        d.seq_size = 0;
        d.set_size = 0;
        for (auto &contrast : d) {
          d.seq_size += contrast.seq_size;
          d.set_size += contrast.set_size;
        }
      }
      break;
    case OccurrenceFilter::MaskOccurrences:
      if (options.verbosity >= Verbosity::verbose)
        cout << "Masking " << motif << " occurrences." << endl;
      for (auto &contrast : d)
        mask_motif_occurrences(motif, contrast, options);
      break;
  }
}
};
