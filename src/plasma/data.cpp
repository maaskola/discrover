/*
 * =====================================================================================
 *
 *       Filename:  data.cpp
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
#include <map>
#include "../aux.hpp"
#include "data.hpp"
#include "../verbosity.hpp"
#include "../sha1.hpp"

using namespace std;

string sha1hash(const string &s)
{
  const char *t = s.c_str();
  unsigned char hash[20];
  char hexstring[41]; // 40 chars + a zero
  int end = (int) strlen(t);
  sha1::calc(t, end, hash);
  sha1::toHexString(hash, hexstring);
  string res;
  for(size_t i = 0; i < 40; i++)
    res += hexstring[i];
  return(res);
}

namespace Seeding {
  DataSet::DataSet() :
  Data::Basic::Set<Fasta::Entry>() {
  }

  DataSet::DataSet(const Specification::DataSet &s, bool revcomp, size_t n_seq) :
    Data::Basic::Set<Fasta::Entry>(s, revcomp, n_seq) { };

  RemovalReport DataSet::mask(const Seeding::mask_sub_t &mask, char mask_symbol) {
    RemovalReport report;
    Verbosity verbosity = Verbosity::info;
    for(auto &seq: sequences) {
      auto iter = mask.find(seq.definition);
      if(iter != end(mask)) {
        report.sequences++;
        if(verbosity >= Verbosity::debug) {
          cout << "Masking " << path << " / " << seq.definition << endl
            << seq.sequence << endl;
          size_t idx = 0;
          for(auto pos: iter->second) {
            while(idx < pos) {
              cout << "-";
              idx++;
            }
            cout << "X";
            idx++;
            report.nucleotides++;
          }
          while(idx < seq.sequence.size()) {
            cout << "-";
            idx++;
          }
          cout << endl;
        }
        for(auto pos: iter->second)
          seq.sequence[pos] = mask_symbol;
        if(verbosity >= Verbosity::debug)
          if(iter != end(mask))
            cout << seq.sequence << endl;
      }
    }
    return(report);
  }

  RemovalReport DataSet::drop_sequences(const unordered_set<string> &ids) {
    const bool noisy_output = false;

    RemovalReport report;
    auto iter = sequences.rbegin();
    size_t idx = 0;
    while(iter != sequences.rend()) {
      if(noisy_output) {
        cerr << "idx = " << idx << endl;
        cerr << "Checking " << path << " " << iter->definition << " for dropping" << endl;
      }
      if(ids.find(iter->definition) != end(ids)) {
        if(noisy_output)
          cerr << "Dropping!" << endl;
        report.sequences++;
        report.nucleotides += iter->sequence.size();
        set_size--;
        seq_size -= iter->sequence.size(); // TODO: find out if this is correct for reverse complements
        auto i = iter;
        bool done = (++i) == sequences.rend();
        sequences.erase(--(iter++).base());
        if(done) break;
      } else
        iter++;
      if(noisy_output)
        cout << "idxB = " << idx++ << endl;
    }
    return(report);
  }

  DataSeries::DataSeries() :
  Data::Basic::Series<DataSet>() {
  }

  DataSeries::DataSeries(const std::string &name, const Specification::DataSets &paths, bool revcomp, size_t n_seq) :
  Data::Basic::Series<DataSet>(name, paths, revcomp, n_seq) {
  }

  DataCollection::DataCollection(const Specification::DataSets &paths, bool revcomp, size_t n_seq) :
  Data::Basic::Collection<DataSeries>(paths, revcomp, n_seq) {
  }

  RemovalReport DataSeries::mask(const mask_t &mask) {
    RemovalReport report;
    for(auto &data_set: sets) {
      auto iter = mask.find(data_set.path);
      if(iter != end(mask))
        report += data_set.mask(iter->second);
    }
    return(report);
  }

  RemovalReport DataSeries::drop_sequences(map<string, unordered_set<string>> &ids) {
    RemovalReport report;
    for(auto &data_set: sets) {
      auto iter = ids.find(data_set.path);
      if(iter != end(ids))
        report += data_set.drop_sequences(iter->second);
    }
    return(report);
  }

  RemovalReport DataCollection::mask(const mask_t &mask) {
    RemovalReport report;
    for(auto &data_series: series)
      report += data_series.mask(mask);
    return(report);
  }

  RemovalReport DataCollection::drop_sequences(map<string, unordered_set<string>> &ids) {
    RemovalReport report;
    for(auto &data_series: series)
      report += data_series.drop_sequences(ids);
    return(report);
  }

  RemovalReport::RemovalReport(size_t n, size_t s) : nucleotides(0), sequences(0) { };
  RemovalReport &operator+=(RemovalReport &a, const RemovalReport &b) {
    a = a + b;
    return(a);
  }
  RemovalReport operator+(const RemovalReport &a, const RemovalReport &b) {
    return(RemovalReport(a.nucleotides + b.nucleotides, a.sequences + b.sequences));
  }
}

