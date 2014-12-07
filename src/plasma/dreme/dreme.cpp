#include "dreme.hpp"
#include "../../hmm/basedefs.hpp"
#include "../../dreme_config.hpp"
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

namespace Dreme {
const char *BinaryNotFoundException::what() const throw() {
  return "Could not find Dreme binary.";
};

InvalidLengthsException::InvalidLengthsException(size_t s1, size_t s2)
    : exception(), min_size(s1), max_size(s2){};
const char *InvalidLengthsException::InvalidLengthsException::what() const
    throw() {
  return (("min_size (" + boost::lexical_cast<string>(min_size)
           + ") is not smaller than max_size ("
           + boost::lexical_cast<string>(max_size) + ").").c_str());
};

list<pair<string, double>> parse_dreme_output(const string &dir) {
  list<pair<string, double>> motifs;

  string path = dir + "/dreme.txt";

  ifstream ifs(path.c_str());

  string line;
  while (getline(ifs, line))
    if (line.substr(0, 6) == "MOTIF ") {
      string motif = line.substr(6);
      double score = 0;  // FIXME
      motifs.push_back(make_pair(motif, score));
    }
  return (motifs);
}

list<pair<string, double>> run(const string &path1, const string &path2,
                               size_t min_size, size_t max_size, bool revcomp,
                               size_t n_motifs, bool remove_temp_dir) {
  if (path2 == "")
    cout << "Running DREME for one FASTA file: " << path1 << endl;
  else
    cout << "Running DREME for two FASTA files: " << path1 << " " << path2
         << endl;

  stringstream str;
  str << DREME_PATH << " -p " << path1;

  if (path2 != "")
    str << " -n " << path2;

  if (min_size == max_size) {
    size_t len = min_size;
    if (len != 0)
      str << " -k " << len;
  } else {
    if (min_size > max_size)
      throw InvalidLengthsException(min_size, max_size);
    str << " -mink " << min_size << " -maxk " << max_size;
  }

  if (not revcomp)
    str << " -norc";

  if (n_motifs > 0)
    str << " -m " << n_motifs;

  auto out_path = boost::filesystem::temp_directory_path()
                  / "dreme-out-%%%%-%%%%";
  ;
  string dreme_output_dir
      = boost::filesystem::unique_path(out_path.string()).string();

  str << " -oc " << dreme_output_dir;

  string command = str.str();

  cout << "Command for running DREME = " << command << endl;

  int res = system(command.c_str());
  if (res != 0) {
    cout << "There was a problem executing DREME. Exiting." << endl;
    exit(res);
  }

  auto regexes = parse_dreme_output(dreme_output_dir);

  if (remove_temp_dir) {
    cout << "Removing temporary directory for DREME output " << dreme_output_dir
         << endl;
    boost::filesystem::remove_all(dreme_output_dir);
  }

  return (regexes);
}
}
