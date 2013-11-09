/* =====================================================================================
 * Copyright (c) 2013, Jonas Maaskola
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 *
 *       Filename:  newanalysis.hpp
 *
 *    Description:  A new analysis to find seeds for motifs
 *
 *        Created:  Sat Nov 09 15:17:21 2013 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <ctime>          // for clock()
#include <sys/times.h>    // for times(struct tms *)
#include <boost/program_options.hpp>
#include "align.hpp"
#include "specification.hpp"
#include "data.hpp"
#include "../GitSHA1.hpp"
#include "../verbosity.hpp"

const std::string header = "# How to interpret this file:\n"
"# The program proceeds iteratively, at each step determining the single most discriminative word.";

const std::string program_name = "newanalysis";

std::string gen_usage_string()
{
  const std::string usage = "This program executes a progressive algorithm, which in each iteration\n"
    "determines the IUPAC word with the highest residual score that is enriched,\n"
    "in signal.fa, accepts it, and masks all occurrences or removes all data samples\n"
    "in which it occurs, before proceeding.\n"
    "\n"
    "Some example calls:\n" +
    program_name + " -f signal.fa -f control.fa -m 8 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 -d 0 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 -s freq -d 0 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 -s freq -d 2 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 -d 2 > analysis.txt\n" +
    program_name + " -f signal.fa -f control.fa -m 4-12 --rdeg 0.2 > analysis.txt\n" +
    program_name + " -f signal:signal.fa -f control.fa -m signal:8 > analysis.txt\n";
  return usage;
}

namespace Analysis {
  struct options_t {
    Verbosity verbosity;
    Specification::DataSets paths;
    size_t max_seed_length;
    bool revcomp;
    size_t n_seq;
  };
  std::ostream &operator<<(std::ostream &os, const options_t& options) {
    os << "Verbosity = " << static_cast<int>(options.verbosity) << std::endl;
    for(auto &path: options.paths)
      os << "Path = " << path << std::endl;
    return(os);
  }
}

bool biggerScore(const NmerStats &a, const NmerStats &b) {
  return(a.score > b.score);
}

std::ostream &operator<<(std::ostream &os, const NmerStats &stats) {
  os << stats.nmer << " " << stats.score;
  return(os);
}

int perform_analysis(const Analysis::options_t &options) {
  if(options.verbosity >= Verbosity::verbose)
    for(auto &path: options.paths)
      std::cout << "Using path: " << path << std::endl;

  Plasma::DataCollection data_collection(options.paths, options.revcomp, options.n_seq);

  NucleotideIndex<size_t, size_t> index(data_collection, options.verbosity);

  std::vector<NmerStats> nmers;
  for(size_t k = 1; k <= options.max_seed_length; k++)
    for(auto &nmer: index.nmer_analysis(k, options.verbosity))
      nmers.push_back(nmer);

  std::sort(begin(nmers), end(nmers), biggerScore);
  for(auto &nmer: nmers)
    std::cout << nmer << std::endl;

  return EXIT_SUCCESS;
}

using namespace std;

int main(int argc, const char** argv) {
  clock_t start_time, end_time;
  double cpu_time_used;
  start_time = clock();

  Analysis::options_t options;

  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Basic options");
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("verbose,v", "Be verbose about the progress")
    ("noisy,V", "Be very verbose about the progress")
    ;

  po::positional_options_description pos;
  pos.add("fasta",-1);

  desc.add_options()
    ("fasta,f", po::value<Specification::DataSets>(&options.paths)->required(), "FASTA file(s) with nucleic acid sequences.")
    ("k", po::value<size_t>(&options.max_seed_length)->default_value(8), "Consider seeds up to this length. Default = 8")
    ("revcomp,r", po::bool_switch(&options.revcomp), "Also consider the reverse complements of the sequences.")
    ("nseq", po::value<size_t>(&options.n_seq)->default_value(0), "Use only the first N sequences of each file. Use 0 to indicate all sequences.")
    ;
 
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
    options.verbosity = Verbosity::verbose;
  if(vm.count("noisy"))
    options.verbosity = Verbosity::debug;

  if(vm.count("version") and not vm.count("help"))
  {
    cout << program_name << " " << GIT_DESCRIPTION << endl;
    if(options.verbosity >= Verbosity::info)
      cout << GIT_SHA1 << endl;
    return(EXIT_SUCCESS);
  }

  if (vm.count("help")) {
    cout << program_name << " " << GIT_DESCRIPTION << endl << "Copyright (C) 2011 Jonas Maaskola\n"
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

  int return_value = perform_analysis(options);

  end_time = clock();
  cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

  if(options.verbosity >= Verbosity::info) {
    cout << "CPU time used = " << cpu_time_used << endl;
    if(options.verbosity >= Verbosity::verbose) {
      tms tm;
      clock_t some_time = times(&tm);
      cout << "times() return " << ((double) some_time) << endl;
      cout << "utime = " << tm.tms_utime << endl;
      cout << "stime = " << tm.tms_stime << endl;
      cout << "cutime = " << tm.tms_cutime << endl;
      cout << "cstime = " << tm.tms_cstime << endl;
    }
  }

  return return_value;
}
