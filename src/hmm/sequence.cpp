#include <fstream>
#include "../plasma/fasta.hpp"
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "sequence.hpp"
#include "../plasma/io.hpp"

using namespace std;

vector<string> extract_seq_ids(const string &path, size_t nseq,
                               Verbosity verbosity) {
  if (verbosity >= Verbosity::info)
    cout << "Getting sequence IDs from " << path << endl;
  if (not boost::filesystem::exists(path))
    throw Exception::File::Existence(path);
  vector<string> s;
  ifstream f(path.c_str());
  size_t idx = 0;
  auto parsing = [&](Fasta::Entry &&entry) {
    s.push_back(move(entry.definition));
    return nseq == 0 or idx++ < nseq;
  };
  auto parser = Fasta::make_parser(parsing);
  f >> parser;
  f.close();
  return s;
}

string seq2string(const seq_t &seq) {
  string s;
  char c;
  for (size_t i = 0; i < seq.size(); i++) {
    switch (seq[i]) {
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
  return s;
}
