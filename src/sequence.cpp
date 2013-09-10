#include <fstream>
#include "plasma/fasta.hpp"
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "sequence.hpp"

using namespace std;

vector<string> extract_seq_ids(const string &path, size_t nseq, Verbosity verbosity)
{
  if(verbosity >= Verbosity::info)
    cout << "Getting sequence IDs from " << path << endl;
  if(not boost::filesystem::exists(path))
  {
    cout << "Error: " << path << " does not exist!" << endl;
    exit(-1);
  }
  vector<string> s;
  ifstream f(path.c_str());
  size_t idx = 0;
  auto parsing = [&](Fasta::Entry &&entry) {
    s.push_back(std::move(entry.definition));
    return(nseq == 0 or idx++ < nseq);
  };
  auto parser = Fasta::make_parser(parsing);
  f >> parser;
  f.close();
  return(s);
}

seq_t string2seq(const string &s, int n_enc)
{
  seq_t seq(s.size());
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
        seq[idx] = empty_symbol;
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

string seq2string(const seq_t &seq)
{
  string s;
  char c;
  for(size_t i = 0; i < seq.size(); i++)
  {
    switch(seq[i]) {
      case 0:
        c = 'a';
        break;
      case 1:
        c = 'c';
        break;
      case 2:
        c = 'g';
        break;
      case 3:
        c = 't';
        break;
      case empty_symbol:
        c = '$';
        break;
      default:
        c = 'n';
        break;
    }
    s.push_back(c);
  }
  return(s);
}

seq_t random_seq(size_t n, size_t alphabet_size)
{
  seq_t s(n);
  for(size_t i = 0; i < n; i++)
    s(i) = (rand() % alphabet_size);
  return(s);
}

