#include "hmm/basedefs.hpp"
#include "dreme_config.hpp"
#include <sstream>
#include <list>
#include <fstream>
#include <exception>

using namespace std;


namespace Dreme {
  class BinaryNotFoundException: public exception {
    virtual const char* what() const throw() {
      return "Could not find Dreme binary.";
    }
  };

  class InvalidLengthsException: public exception {
    size_t min_size, max_size;
    public:
    InvalidLengthsException(size_t s1, size_t s2) : exception(), min_size(s1), max_size(s2) { };
    virtual const char* what() const throw()
    {
      return(("min_size (" + boost::lexical_cast<string>(min_size) + ") is not smaller than max_size (" + boost::lexical_cast<string>(max_size) + ").").c_str());
    }
  };

  list<string> parse_dreme_output(const string &dir) {
    list<string> motifs;

    string path = dir + "/dreme.txt";

    ifstream ifs(path.c_str());

    string line;
    while(getline(ifs, line))
      if(line.substr(0, 6) == "MOTIF ")
        motifs.push_back(line.substr(6));
    return(motifs);
  }

  list<string> run(const string &path1,
      const string &path2="",
      size_t min_size=0,
      size_t max_size=0,
      const string &path="") {
    if(path2 == "")
      cout << "Running DREME for one FASTA file: " << path1 << endl;
    else
      cout << "Running DREME for two FASTA files: " << path1 << " " << path2 << endl;

    stringstream str;
    str << DREME_PATH << " -p " << path1;

    if(path2 != "")
      str << " -n " << path2;

    if(min_size == max_size) {
      size_t len = min_size;
      if(len != 0)
        str << " -k " << len;
    } else {
      if(min_size > max_size)
        throw InvalidLengthsException(min_size, max_size);
      str << " -mink " << min_size << " -maxk " << max_size;
    }

    string dreme_output_dir;
    if(path == "")
      dreme_output_dir = "/tmp/dlhmm_dreme_out/";
    else
      dreme_output_dir = "/tmp/" + path;

    str << " -oc " << dreme_output_dir;

    string command = str.str();

    cout << "Command = " << command << endl;

    int res = system(command.c_str());

    list<string> regexes = parse_dreme_output(dreme_output_dir);
    for(auto &x: regexes)
      cout << "Dreme found motif: " << x << endl;
    return(regexes);
  }
}

int main(int argc, const char** argv) {
  if(argc < 2) {
    cout << "Please provide at least one path to a FASTA file." << endl;
    return -1;
  }
  if(argc > 3) {
    cout << "Please provide at most two paths to a FASTA file." << endl;
    return -2;
  }

  string path1 = argv[1];
  string path2 = "";

  if(argc == 3)
    path2 = argv[2];

  /*
  try {
    auto regexes = Dreme::run(path1, path2);
  } catch (exception& e) {
    cout << "Exception: " << e.what() << endl;
  } */
  try {
    auto regexes = Dreme::run(path1, path2, 8, 8);
  } catch (exception& e) {
    cout << "Exception: " << e.what() << endl;
  }
  /*
  try {
    auto regexes = Dreme::run(path1, path2, 6, 8);
  } catch (exception& e) {
    cout << "Exception: " << e.what() << endl;
  }
  try {
    auto regexes = Dreme::run(path1, path2, 8, 6);
  } catch (exception& e) {
    cout << "Exception: " << e.what() << endl;
  }  */

  return EXIT_SUCCESS;
}
