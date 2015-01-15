
#include <random>
#include "fasta.hpp"
#include "../shuffle/dinucleotide_shuffle.hpp"
#include "../aux.hpp"
#include "io.hpp"
#include "../verbosity.hpp"

using namespace std;

namespace Fasta {

/* Valid nucleic acid codes according to
 * http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
 *  A  adenosine          C  cytidine             G  guanine
 *  T  thymidine          N  A/G/C/T (any)        U  uridine
 *  K  G/T (keto)         S  G/C (strong)         Y  T/C (pyrimidine)
 *  M  A/C (amino)        W  A/T (weak)           R  G/A (purine)
 *  B  G/T/C              D  G/A/T                H  A/C/T
 *  V  G/C/A              -  gap of indeterminate length
 */
// Disallow gaps of indeterminate length '-'
static const string valid_nucleotides = "acgtnuksymwrbdhy";

mt19937 Fasta::EntropySource::shuffling_rng;
mt19937 Fasta::EntropySource::random_nucl_rng;
uniform_int_distribution<size_t> r_unif;
uniform_int_distribution<size_t> r_nucl(0, 3);

Entry::Entry() : definition(), sequence(){};
Entry::Entry(const Entry &entry)
    : definition(entry.definition), sequence(entry.sequence){};
Entry::Entry(const IEntry &ientry) : definition(ientry.definition), sequence() {
  size_t pos = ientry.sequence.find("$");
  sequence = ientry.sequence.substr(0, pos);
}

size_t Entry::mask(const vector<size_t> &positions) {
  size_t masked_nucleotides = 0;
  const char mask_symbol = 'n';
  const Verbosity verbosity = Verbosity::info;
  if (verbosity >= Verbosity::debug) {
    cout << "Masking " << definition << endl << sequence << endl;
    size_t idx = 0;
    for (auto pos : positions) {
      while (idx < pos) {
        cout << "-";
        idx++;
      }
      cout << "X";
      idx++;
      masked_nucleotides++;
    }
    while (idx < sequence.size()) {
      cout << "-";
      idx++;
    }
    cout << endl;
  }
  for (auto pos : positions) {
    if (pos >= sequence.size())
      pos = 2 * sequence.size() - pos;
    sequence[pos] = mask_symbol;
  }
  if (verbosity >= Verbosity::debug)
    cout << "Masked  " << definition << endl << sequence << endl;
  return masked_nucleotides;
}

/** Convert a string of nucleic acids into an IEntry::seq_t
 * Note: characters other than 'acgt' are assigned random nucleotides
 *       if n_enc < 0.
 **/
IEntry::seq_t string2seq(const string &s, int n_enc = -1) {
  IEntry::seq_t seq(s.size());
  size_t idx = 0;
  for (auto iter : s) {
    switch (tolower(iter)) {
      case 'a':
        seq[idx] = 0;
        break;
      case 'c':
        seq[idx] = 1;
        break;
      case 'g':
        seq[idx] = 2;
        break;
      case 't':
      case 'u':
        seq[idx] = 3;
        break;
      case '$':
        seq[idx] = IEntry::empty_symbol;
        break;
      default:
        if (n_enc < 0)
          seq[idx] = r_nucl(EntropySource::random_nucl_rng);
        else
          seq[idx] = n_enc;
    }
    idx++;
  }
  return seq;
}

IEntry::IEntry(const Entry &entry)
    : Entry(entry), isequence(string2seq(entry.sequence)){};

size_t IEntry::mask(const vector<size_t> &positions) {
  size_t masked_nucleotides = 0;
  // uses random nucleotides;
  // const char mask_symbol = 'n';
  // const Verbosity verbosity = Verbosity::info;
  size_t dollar_position = sequence.find("$");
  for (auto pos : positions) {
    // TODO do something more sensible than random nucleotides
    size_t nucl_idx = r_nucl(EntropySource::random_nucl_rng);
    sequence[pos] = "acgt"[nucl_idx];
    if (dollar_position != string::npos)
      sequence[sequence.size() - pos - 1] = "tgca"[nucl_idx];
  }

  isequence = string2seq(sequence);
  return masked_nucleotides;
}

ostream &operator<<(ostream &os, const Entry &entry) {
  os << entry.string();
  return os;
}

istream &operator>>(istream &is, Entry &entry) {
  // consume until '>'
  char c;
  while (is >> c)
    if (c == '>')
      break;

  if (not is.good())
    return is;

  // read the definition
  getline(is, entry.definition);

  if (is.eof() or not is.good()) {
    is.setstate(ios_base::failbit);
    return is;
  }

  // consume the sequence
  entry.sequence = "";
  while (is >> c) {
    if (c == '>') {
      is.putback(c);
      break;
    } else if (valid_nucleotides.find(tolower(c)) != string::npos)
      entry.sequence += c;
    else if (not isspace(c))
      throw Exception::NucleicAcids::InvalidNucleotideCode(c);
  }

  // the while loop may end because of reaching end-of-file
  // in this case we want to clear the fail bit
  is.clear((is.eof() ? ios_base::eofbit : ios_base::goodbit)
         | (is.bad() ? ios_base::badbit : ios_base::goodbit));

  return is;
}

istream &operator>>(istream &is, IEntry &ientry) {
  Entry entry;
  if (is >> entry)
    ientry = IEntry(entry);
  return is;
}

istream &operator>>(istream &is, vector<Entry> &seqs) {
  Entry entry;
  while (is >> entry)
    seqs.push_back(move(entry));

  // vector<Entry> is defined as zero or more sequences
  // so if no entries are read and end-of-file is reached, we clear the fail bit
  is.clear((is.eof() ? ios_base::eofbit : ios_base::goodbit)
         | (is.bad() ? ios_base::badbit : ios_base::goodbit));

  return is;
}

istream &operator>>(istream &is, vector<IEntry> &seqs) {
  IEntry entry;
  while (is >> entry)
    seqs.push_back(move(entry));

  // vector<IEntry> is defined as zero or more sequences
  // so if no entries are read and end-of-file is reached, we clear the fail bit
  is.clear((is.eof() ? ios_base::eofbit : ios_base::goodbit)
         | (is.bad() ? ios_base::badbit : ios_base::goodbit));

  return is;
}

void read_fasta(const string &path, vector<Entry> &sequences, bool revcomp,
                size_t n_seq, bool shuffled) {
  try {
    parse_file(path, [&](istream &is) { is >> sequences; });
  } catch (runtime_error &e) {
    std::cout << "Error while reading FASTA file " << path << "." << std::endl;
    throw e;
  }

  if (n_seq > 0 and sequences.size() > n_seq)  // TODO: improve efficiency by
                                               // only reading in n_seq
                                               // sequences
    sequences.resize(n_seq);

  if (sequences.size() == 0)
    // TODO: throw exception?
    cout << "Warning: no sequences found while parsing " << path << " as FASTA format file." << endl
         << "Please check the format of this file or whether it is the correct one." << endl;

  if (shuffled)
    for (auto &s : sequences) {
      s.definition = "Shuffle of " + s.definition;
      s.sequence = dinucleotideShuffle(s.sequence,
                                       r_unif(EntropySource::shuffling_rng));
    }
};

void read_fasta(const string &path, vector<IEntry> &isequences, bool revcomp,
                size_t n_seq, bool shuffled) {
  vector<Entry> sequences;
  read_fasta(path, sequences, revcomp, n_seq, shuffled);
  for (auto &s : sequences) {
    IEntry is(s);
    if (revcomp)
      is.sequence += "$" + reverse_complement(is.sequence);
    is.isequence = string2seq(is.sequence);
    isequences.push_back(is);
  }
};
}

string reverse_complement(const string &s) {
  string t;
  for (string::const_reverse_iterator iter = s.rbegin(); iter != s.rend();
       iter++)
    switch (*iter) {
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
        throw Exception::NucleicAcids::InvalidNucleotideCode(*iter);
    }
  return t;
}
