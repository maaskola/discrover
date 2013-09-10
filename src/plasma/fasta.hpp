
#ifndef FASTA_HPP
#define FASTA_HPP

#include <iostream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

std::string reverse_complement(const std::string &s);

namespace Fasta {
  struct Entry {
    std::string definition;
    std::string sequence;
    std::string string(size_t width=60) const { return(">" + definition + "\n" + sequence); };
    std::string reverse_complement() const { return(::reverse_complement(sequence)); };
  };
  struct IEntry : public Entry {
    typedef unsigned char alphabet_idx_t;
    typedef boost::numeric::ublas::vector<alphabet_idx_t> seq_t;
    // typedef std::vector<alphabet_idx_t> seq_t;
    static const alphabet_idx_t empty_symbol = 5;
    seq_t isequence;
    IEntry(const Entry &entry=Entry());
  };

  template <typename X>
    struct Parser {
      X x;
      bool operator()(Entry &&entry) {
        return(x(std::move(entry)));
      };
    };
  template <typename T>
    Parser<typename std::decay<T>::type>
    make_parser(T&& t) {
      return { std::forward<T>(t) };
    }
};

std::istream &operator>>(std::istream &is, Fasta::Entry &entry);
std::istream &operator>>(std::istream &is, Fasta::IEntry &entry);

template <class X>
std::istream &operator>>(std::istream &is, Fasta::Parser<X> &parser) {
  bool proceed = true;
  while(proceed and is.good()) {
    Fasta::Entry entry;
    is >> entry;
    proceed = parser(std::move(entry));
  }
  return(is);
};

std::ostream &operator<<(std::ostream &os, const Fasta::Entry &entry);
std::istream &operator>>(std::istream &is, std::vector<Fasta::Entry> &parser);
std::istream &operator>>(std::istream &is, std::vector<Fasta::IEntry> &parser);

void read_fasta(const std::string &path, std::vector<Fasta::Entry> &entries, bool revcomp, size_t n_seq=0);
void read_fasta(const std::string &path, std::vector<Fasta::IEntry> &entries, bool revcomp, size_t n_seq=0);

#endif

