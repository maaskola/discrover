
#include "fasta.hpp"
#include "io.hpp"

using namespace std;

Fasta::IEntry::seq_t string2seq_(const string &s, int n_enc=-1)
{
  Fasta::IEntry::seq_t seq(s.size());
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
        seq[idx] = Fasta::IEntry::empty_symbol;
      default:
        if(n_enc < 0)
          seq[idx] = random() % 4;
        else
          seq[idx] = n_enc;
        // throw("Wrong encoding");
    }
    idx++;
  }
  return(seq);
}

namespace Fasta {
  IEntry::IEntry(const Entry &entry) : Entry(entry), isequence(string2seq_(entry.sequence)) {
  };
}

ostream &operator<<(ostream &os, const Fasta::Entry &entry) {
  os << entry.string();
  return(os);
}

istream &operator>>(istream &is, Fasta::Entry &entry) {
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

istream &operator>>(istream &is, Fasta::IEntry &ientry) {
  Fasta::Entry entry;
  is >> entry;
  ientry = Fasta::IEntry(entry);
  return(is);
}

istream &operator>>(istream &is, vector<Fasta::Entry> &seqs) {
  while(is.good()) {
    Fasta::Entry entry;
    is >> entry;
    seqs.push_back(move(entry));
  }
  return(is);
}

istream &operator>>(istream &is, vector<Fasta::IEntry> &seqs) {
  while(is.good()) {
    Fasta::IEntry entry;
    is >> entry;
    seqs.push_back(move(entry));
  }
  return(is);
}

void read_fasta(const string &path, vector<Fasta::Entry> &sequences, bool revcomp, size_t n_seq) {
  parse_file(path, [&](istream &is) { is >> sequences; });
  if(n_seq > 0) // TODO: improve efficiency by only reading in n_seq sequences
    sequences.resize(n_seq);
};

void read_fasta(const string &path, vector<Fasta::IEntry> &sequences, bool revcomp, size_t n_seq) {
  parse_file(path, [&](istream &is) { is >> sequences; });
  if(n_seq > 0) // TODO: improve efficiency by only reading in n_seq sequences
    sequences.resize(n_seq);
  if(revcomp)
    for(auto &s: sequences) {
      s.sequence += "$" + reverse_complement(s.sequence);
      s.isequence = string2seq_(s.sequence);
    }
};

using namespace std;

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

