
#include <random>
#include "fasta.hpp"
#include "../shuffle/dinucleotide_shuffle.hpp"
#include "io.hpp"
#include "../verbosity.hpp"

using namespace std;

namespace Fasta {
  mt19937 Fasta::EntropySource::shuffling_rng;
  mt19937 Fasta::EntropySource::random_nucl_rng;
  uniform_int_distribution<size_t> r_unif;

  Entry::Entry() : definition(), sequence() { };
  Entry::Entry(const Entry &entry) : definition(entry.definition), sequence(entry.sequence) { };
  Entry::Entry(const IEntry &ientry) : definition(ientry.definition), sequence()
  {
    size_t pos = ientry.sequence.find("$");
    sequence = ientry.sequence.substr(0,pos);
  }

  size_t Entry::mask(vector<size_t> positions) {
    size_t masked_nucleotides = 0;
    const char mask_symbol = 'n';
    const Verbosity verbosity = Verbosity::info;
    if(verbosity >= Verbosity::debug) {
      cout << "Masking " << definition << endl << sequence << endl;
      size_t idx = 0;
      for(auto pos: positions) {
        while(idx < pos) {
          cout << "-";
          idx++;
        }
        cout << "X";
        idx++;
        masked_nucleotides++;
      }
      while(idx < sequence.size()) {
        cout << "-";
        idx++;
      }
      cout << endl;
    }
    for(auto pos: positions) {
      if(pos >= sequence.size())
        pos = 2 * sequence.size() - pos;
      sequence[pos] = mask_symbol;
    }
    if(verbosity >= Verbosity::debug)
      cout << "Masked  " << definition << endl << sequence << endl;
    return(masked_nucleotides);
  }

  IEntry::seq_t string2seq(const string &s, int n_enc=-1)
  {
    std::uniform_int_distribution<size_t> r_nucl(0,3);
    IEntry::seq_t seq(s.size());
    size_t idx = 0;
    for(auto iter : s)
    {
      switch(iter) {
        case 'a':
        case 'A':
          seq[idx] = 0;
          break;
        case 'c':
        case 'C':
          seq[idx] = 1;
          break;
        case 'g':
        case 'G':
          seq[idx] = 2;
          break;
        case 't':
        case 'T':
        case 'u':
        case 'U':
          seq[idx] = 3;
          break;
        case '$':
          seq[idx] = IEntry::empty_symbol;
          break;
        default:
          if(n_enc < 0)
            seq[idx] = r_nucl(EntropySource::random_nucl_rng);
          else
            seq[idx] = n_enc;
          // throw("Wrong encoding");
      }
      idx++;
    }
    return(seq);
  }

  IEntry::IEntry(const Entry &entry) : Entry(entry), isequence(string2seq(entry.sequence)) {
  };

  size_t IEntry::mask(vector<size_t> positions) {
    size_t masked_nucleotides = 0;
    // uses random nucleotides;
    // const char mask_symbol = 'n';
    // const Verbosity verbosity = Verbosity::info;
    const Verbosity verbosity = Verbosity::debug;
    if(verbosity >= Verbosity::debug) {
      cout << "Masking IEntry " << definition << endl << sequence << endl;
      size_t idx = 0;
      for(auto pos: positions) {
        while(idx < pos) {
          cout << "-";
          idx++;
        }
        cout << "X";
        idx++;
        masked_nucleotides++;
      }
      while(idx < sequence.size()) {
        cout << "-";
        idx++;
      }
      cout << endl;
    }
    size_t dollar_position = sequence.find("$");
    // cout << "dollar_position = " << dollar_position << endl;
    // cout << "string::npos = " << string::npos << endl;
    for(auto pos: positions) {
      size_t nucl_idx = rand() % 4; // TODO do something more sensible than random nucleotides
      sequence[pos] = "acgt"[nucl_idx];
      // cout << "Masked pos " << pos << endl;
      if(dollar_position != string::npos) {
        // if(pos > dollar_position)
        //   sequence[sequence.size() - pos] = "tgca"[nucl_idx];
        // else
        sequence[sequence.size() - pos - 1] = "tgca"[nucl_idx];
        // cout << "Masked pos " << (sequence.size() - pos - 1) << endl;
      }
      // if(pos >= sequence.size())
      //   pos = 2 * sequence.size() - pos;
      // sequence[pos] = mask_symbol;
      // sequence[pos] = "acgt"[rand() % 4]; // TODO do something more sensible than random nucleotides
    }

    isequence = string2seq(sequence);
    if(verbosity >= Verbosity::debug)
      cout << sequence << endl;
      // cout << "Masked  " << definition << endl << sequence << endl;
    return(masked_nucleotides);
  }



  ostream &operator<<(ostream &os, const Entry &entry) {
    os << entry.string();
    return(os);
  }

  istream &operator>>(istream &is, Entry &entry) {
    if(not is.good())
      return(is);
    // consume until '>'
    char c;
    is.get(c);
    while(c != '>') {
      if(not is.good())
        return(is);
      is.get(c);
    }

    // read the definition
    getline(is, entry.definition);

    // consume the sequence
    entry.sequence = "";
    while(is.good()) {
      is.get(c);
      if(c == '>') {
        is.putback(c);
        break;
      } else if(c >= ' ' and c <= '~') // the least and largest printable ASCII codes
        entry.sequence += c;
    }
    return(is);
  }

  istream &operator>>(istream &is, IEntry &ientry) {
    Entry entry;
    is >> entry;
    ientry = IEntry(entry);
    return(is);
  }

  istream &operator>>(istream &is, vector<Entry> &seqs) {
    while(is.good()) {
      Entry entry;
      is >> entry;
      seqs.push_back(move(entry));
    }
    return(is);
  }

  istream &operator>>(istream &is, vector<IEntry> &seqs) {
    while(is.good()) {
      IEntry entry;
      is >> entry;
      seqs.push_back(move(entry));
    }
    return(is);
  }

  void read_fasta(const string &path, vector<Entry> &sequences, bool revcomp, size_t n_seq, bool shuffled) {
    parse_file(path, [&](istream &is) { is >> sequences; });
    if(n_seq > 0 and sequences.size() > n_seq) // TODO: improve efficiency by only reading in n_seq sequences
      sequences.resize(n_seq);
    if(shuffled)
      for(auto &s: sequences)
        s.sequence = dinucleotideShuffle(s.sequence, r_unif(EntropySource::shuffling_rng));
  };

  void read_fasta(const string &path, vector<IEntry> &isequences, bool revcomp, size_t n_seq, bool shuffled) {
    vector<Entry> sequences;
    read_fasta(path, sequences, revcomp, n_seq, shuffled);
    for(auto &s: sequences) {
      IEntry is(s);
      if(revcomp)
        is.sequence += "$" + reverse_complement(is.sequence);
        is.isequence = string2seq(is.sequence);
      isequences.push_back(is);
    }
  };
}

string reverse_complement(const string &s) {
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

