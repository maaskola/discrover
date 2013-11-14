/*
 * =====================================================================================
 *
 *       Filename:  align.cpp
 *
 *    Description:  Alignment routines based on suffix arrays
 *
 *        Version:  1.0
 *        Created:  14.08.2012 22:23:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola
 *   Organization:  
 *
 * =====================================================================================
 */

#include <limits>
#include "align.hpp"
#include "fasta.hpp"

using namespace std;

bool iupac_included(char r, char q)
{
  q = tolower(q);
  r = tolower(r);
  if(q==r)
    return(true);
  switch(r) {
    case 'a':
      return(q=='w' or q=='m' or q=='r' or q=='d' or q=='h' or q=='v' or q=='n');
    case 'c':
      return(q=='s' or q=='m' or q=='y' or q=='b' or q=='h' or q=='v' or q=='n');
    case 'g':
      return(q=='s' or q=='k' or q=='r' or q=='b' or q=='d' or q=='v' or q=='n');
    case 't':
      return(q=='w' or q=='k' or q=='y' or q=='b' or q=='d' or q=='h' or q=='n');
    default:
      return(false);
  }
}

string iupac_reverse_complement(const string &s) {
  string t;
  for(string::const_reverse_iterator iter = s.rbegin(); iter != s.rend(); iter++)
    switch(*iter) {
      case 'a':
        t += "t";
        break;
      case 'c':
        t += "g";
        break;
      case 'g':
        t += "c";
        break;
      case 't':
      case 'u':
        t += "a";
        break;
      case 'n':
        t += "n";
        break;
      case 'A':
        t += "T";
        break;
      case 'C':
        t += "G";
        break;
      case 'G':
        t += "C";
        break;
      case 'T':
      case 'U':
        t += "A";
        break;
      case 'N':
        t += "N";
        break;
      case 'b':
        t += 'v';
        break;
      case 'v':
        t += 'b';
        break;
      case 'd':
        t += 'h';
        break;
      case 'h':
        t += 'd';
        break;
      case 's':
        t += 's';
        break;
      case 'w':
        t += 'w';
        break;
      case 'm':
        t += 'k';
        break;
      case 'k':
        t += 'm';
        break;
      case 'r':
        t += 'y';
        break;
      case 'y':
        t += 'r';
        break;
      case 'B':
        t += 'V';
        break;
      case 'V':
        t += 'B';
        break;
      case 'D':
        t += 'H';
        break;
      case 'H':
        t += 'D';
        break;
      case 'S':
        t += 'S';
        break;
      case 'W':
        t += 'W';
        break;
      case 'M':
        t += 'K';
        break;
      case 'K':
        t += 'M';
        break;
      case 'R':
        t += 'Y';
        break;
      case 'Y':
        t += 'R';
        break;
      default:
        throw("We have a problem reversing that string!");
    }
  return(t);
}

string read_fasta_with_boundaries(const vector<string> &paths, vector<seq_index_t> &pos2seq, vector<set_index_t> &seq2set, size_t n_seq) {
  string s = "";
  set_index_t set_idx = 0;
  for(auto &path: paths) {
    if(not boost::filesystem::exists(path)) {
      cout << "Error opening file: " << path << " does not exist." << endl;
      exit(-1);
    }

    bool use_gzip = path.substr(path.size() - 3, 3) == ".gz";
    bool use_bzip2 = path.substr(path.size() - 4, 4) == ".bz2";
    ios_base::openmode flags = ios_base::in;
    if(use_gzip or use_bzip2)
      flags |= ios_base::binary;

    ifstream file(path, flags);
    boost::iostreams::filtering_stream<boost::iostreams::input> in;
    if(use_gzip)
      in.push(boost::iostreams::gzip_decompressor());
    if(use_bzip2)
      in.push(boost::iostreams::bzip2_decompressor());
    in.push(file);

    seq_index_t seq_idx = 0;
    auto parsing = [&](Fasta::Entry &&entry) {
      boost::algorithm::to_lower(entry.sequence);
      s += entry.sequence + "$";
      for(size_t i = 0; i < entry.sequence.size() + 1; i++)
        pos2seq.push_back(seq_idx);
      seq2set.push_back(set_idx);
      seq_idx++;
      if(seq_idx == numeric_limits<seq_index_t>::max()) {
        cout << "Sequence index type not large enough." << endl
          << "Please change seq_index_t to a larger type and recompile.";
        throw;
      }
      return(n_seq == 0 or seq_idx < n_seq);
    };
    auto parser = Fasta::make_parser(parsing);
    in >> parser;
    set_idx++;
    if(set_idx == numeric_limits<set_index_t>::max()) {
      cout << "Set index type not large enough." << endl
        << "Please change set_index_t to a larger type and recompile.";
      throw;
    }
  }
  return(s);
}

string collapse_data_series(const Seeding::DataSeries &data_series, vector<seq_index_t> &pos2seq, vector<set_index_t> &seq2set) {
  string s;
  seq_index_t seq_idx = 0;
  set_index_t set_idx = 0;
  for(auto &data_set: data_series) {
    for(auto &seq: data_set) {
      s += seq.sequence + "$";
      for(size_t i = 0; i < seq.sequence.size() + 1; i++)
        pos2seq.push_back(seq_idx);
      seq2set.push_back(set_idx);
      seq_idx++;
      if(seq_idx == numeric_limits<seq_index_t>::max()) {
        cout << "Sequence index type not large enough." << endl
          << "Please change seq_index_t to a larger type and recompile.";
        throw;
      }
    }
    set_idx++;
    if(set_idx == numeric_limits<set_index_t>::max()) {
      cout << "Set index type not large enough." << endl
        << "Please change set_index_t to a larger type and recompile.";
      throw;
    }
  }
  return(s);
}

string collapse_data_collection(const Seeding::DataCollection &collection, vector<seq_index_t> &pos2seq, vector<set_index_t> &seq2set, vector<set_index_t> &set2series) {
  string s;
  seq_index_t seq_idx = 0;
  set_index_t set_idx = 0;
  set_index_t series_idx = 0;
  for(auto &series: collection) {
    for(auto &set: series) {
      for(auto &seq: set) {
        s += seq.sequence + "$";
        for(size_t i = 0; i < seq.sequence.size() + 1; i++)
          pos2seq.push_back(seq_idx);
        seq2set.push_back(set_idx);
        seq_idx++;
        if(seq_idx == numeric_limits<seq_index_t>::max()) {
          cout << "Sequence index type not large enough." << endl
            << "Please change seq_index_t to a larger type and recompile.";
          throw;
        }
      }
      set2series.push_back(series_idx);
      set_idx++;
      if(set_idx == numeric_limits<set_index_t>::max()) {
        cout << "Set index type not large enough." << endl
          << "Please change set_index_t to a larger type and recompile.";
        throw;
      }
    }
    series_idx++;
    if(series_idx == numeric_limits<set_index_t>::max()) {
      cout << "Series index type not large enough." << endl
        << "Please change set_index_t to a larger type and recompile.";
      throw;
    }
  }
  return(s);
}

