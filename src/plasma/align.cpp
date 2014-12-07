/*
 * =====================================================================================
 *
 *       Filename:  align.cpp
 *
 *    Description:  Alignment routines based on suffix arrays
 *
 *        Created:  14.08.2012 22:23:44
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "align.hpp"

using namespace std;

const char TERMINATOR_SYMBOL = '$';

vector<symbol_t> collapse_collection(const Seeding::Collection &collection,
                                     vector<size_t> &pos2seq,
                                     vector<size_t> &seq2set,
                                     vector<size_t> &set2contrast,
                                     bool allow_iupac_wildcards) {
  vector<symbol_t> s;
  size_t seq_idx = 0;
  size_t set_idx = 0;
  size_t contrast_idx = 0;
  for (auto &contrast : collection) {
    for (auto &dataset : contrast) {
      for (auto &seq : dataset) {
        add_sequence(s, seq.sequence + TERMINATOR_SYMBOL,
                     allow_iupac_wildcards);
        for (size_t i = 0; i < seq.sequence.size() + 1; i++)
          pos2seq.push_back(seq_idx);
        seq2set.push_back(set_idx);
        seq_idx++;
      }
      set2contrast.push_back(contrast_idx);
      set_idx++;
    }
    contrast_idx++;
  }
  return s;
}
