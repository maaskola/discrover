/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  A tool to generate shuffles of sequences
 *
 *        Version:  1.0
 *        Created:  19.11.2013 20:54:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */
#include <cstdlib>
#include <fstream>
#include <cstdlib>
#include <cstddef>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>
#include "../GitSHA1.hpp"
#include "../verbosity.hpp"
#include "../plasma/fasta.hpp"
#include "dinucleotide_shuffle.hpp"

const std::string program_name = "dinucleotide_shuffle";

std::string gen_usage_string()
{
  const std::string usage = "Generates dinucleotide frequency preserving shuffles of FASTA files.\n"
    "The code is based on altschulEriksonDinuclShuffle.py by P. Clote, from Oct 2003.\n";
  return usage;
}


using namespace std;

void shuffle(istream &is, size_t n, size_t seed) {
  mt19937 rng;
  rng.seed(seed);
  uniform_int_distribution<size_t> r_unif;
  auto parsing = [&n, &seed, &rng, &r_unif](Fasta::Entry &&entry) {
    for(size_t i = 0; i < n; i++) {
      string seq = entry.sequence;
      for(auto &s: seq) {
        s = tolower(s);
        if(s == 'u')
          s = 't';
      }
      cout << ">" << entry.definition << endl
        << dinucleotideShuffle(seq, r_unif(rng)) << endl;
    }
    return(true);
  };
  auto parser = Fasta::make_parser(parsing);
  is >> parser;
}

int main(int argc, const char **argv)
{
  vector<string> paths;
  size_t n = 1;
  size_t seed = 1;
  Verbosity verbosity = Verbosity::info;

  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("fasta,f", po::value<vector<string>>(&paths), "FASTA file(s) with nucleic acid sequences. If no files are specified, sequences are read from standard input.")
    ("number,n", po::value<size_t>(&n)->default_value(1), "How many shuffles to generate per sequence.")
    ("seed,s", po::value<size_t>(&seed), "Seed to initialize random number generator.")
    ("verbose,v", "Be verbose about the progress")
    ;
 
  po::positional_options_description pos;
  pos.add("fasta",-1);

  po::variables_map vm;

  try {
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  } catch(po::unknown_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " not known." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::ambiguous_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " is ambiguous." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " was specified multiple times." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::multiple_occurrences &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option --" << e.get_option_name() << " was specified multiple times." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
      << "The value specified for option " << e.get_option_name() << " has an invalid format." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Too many positional options were specified." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
      << "Invalid command line syntax." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
      << "There is a programming error related to command line style." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
      << "The configuration file can not be read." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Validation of option " << e.get_option_name() << " failed." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::error &e) {
    cout << "Error while parsing command line options:" << endl
      << "No further information as to the nature of this error is available, please check your command line arguments." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  }

  if(vm.count("verbose"))
    verbosity = Verbosity::verbose;

  if(vm.count("version") and not vm.count("help"))
  {
    cout << program_name << " " << GIT_DESCRIPTION << " [" << GIT_BRANCH << " branch]" << endl;
    if(verbosity >= Verbosity::verbose)
      cout << GIT_SHA1 << endl;
    return(EXIT_SUCCESS);
  }

  if (vm.count("help")) {
    cout << program_name << " " << GIT_DESCRIPTION << endl << "Copyright (C) 2013 Jonas Maaskola\n"
      "Provided under GNU General Public License Version 3 or later.\n"
      "See the file COPYING provided with this software for details of the license.\n" << endl;
    cout << gen_usage_string() << endl;
    cout << desc << "\n";
    return 1;
  }


  try {
    po::notify(vm);
  } catch(po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " was specified multiple times." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
      << "The value specified for option " << e.get_option_name() << " has an invalid format." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Too many positional options were specified." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
      << "Invalid command line syntax." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
      << "There is a programming error related to command line style." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
      << "The configuration file can not be read." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::required_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "The required option " << e.get_option_name() << " was not specified." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Validation of option " << e.get_option_name() << " failed." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  } catch(po::error &e) {
    cout << "Error while parsing command line options:" << endl
      << "No further information as to the nature of this error is available, please check your command line arguments." << endl
      << "Please inspect the command line help with -h or --help." << endl;
    return(-1);
  }

  if(not vm.count("seed"))
    seed = random_device()();


  if(paths.empty()) {
    shuffle(cin, n, seed);
  } else
    for(auto &path: paths) {
      if(boost::filesystem::exists(path)) {
        ifstream ifs(path.c_str());
        shuffle(ifs, n, seed++);
      }
      else {
        cerr << "Error: " << path << " does not exist." << endl;
        exit(-1);
      }
    }

  return EXIT_SUCCESS;
}

