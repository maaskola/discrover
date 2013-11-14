/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:
 *
 * =====================================================================================
 */

#include <ctime>          // for clock()
#include <sys/times.h>    // for times(struct tms *)
#include <iostream>
#include <omp.h>
#include <fstream>
#include <boost/program_options.hpp>
#include "code.hpp"
#include "score.hpp"
#include "align.hpp"
#include "count.hpp"
#include "mask.hpp"
#include "find.hpp"
#include "../timer.hpp"
#include "plasma_cli.hpp"
#include "../GitSHA1.hpp"

const std::string header = "# How to interpret this file:\n"
"# The program proceeds iteratively, at each step determining the single most discriminative word.";

const std::string program_name = "plasma";

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

using namespace std;

int main(int argc, const char** argv) {
  clock_t start_time, end_time;
  double cpu_time_used;
  start_time = clock();

  Seeding::options_t options;

  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Basic options");
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("verbose,v", "Be verbose about the progress")
    ("noisy,V", "Be very verbose about the progress")
    ;
  po::options_description ext_options = gen_iupac_options_description(options);

  desc.add(ext_options);

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
    options.verbosity = Verbosity::verbose;
  if(vm.count("noisy"))
    options.verbosity = Verbosity::debug;

  if(vm.count("version") and not vm.count("help"))
  {
    cout << program_name << " " << GIT_DESCRIPTION << " [" << GIT_BRANCH << " branch]" << endl;
    if(options.verbosity >= Verbosity::verbose)
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


  if(options.motif_specifications.size() == 0) {
    cout << "Error: you must specify motif lengths of interest with the -m switch." << endl;
    exit(-1);
  }

  if(vm.count("threads"))
    omp_set_num_threads(options.n_threads);


  if(options.verbosity >= Verbosity::verbose) {
    cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
    cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
    cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
  }

  Specification::harmonize(options.motif_specifications, options.paths, options.objectives);

  if(options.verbosity >= Verbosity::verbose) {
    cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
    cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
    cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
  }


  Seeding::Plasma plasma(options);
  Seeding::DataCollection ds = plasma.collection;
  typedef Seeding::Result res_t;
  vector<res_t> results;

  size_t n = options.motif_specifications.size();
  for(auto &motif: options.motif_specifications) {
    for(auto &r:  plasma.find(motif, options.objectives, false))
      results.push_back(r);
    if(--n > 0)
      plasma.apply_mask(results);
  }

  if(results.size() == 1)
    report(cout, results[0], ds, options);
  else {
    sort(begin(results), end(results), [](const res_t &a, const res_t &b) { return(a.log_p <= b.log_p); });
    Seeding::DataCollection original_ds = ds;
    Seeding::options_t opts = options;
    opts.occurrence_filter = Seeding::OccurrenceFilter::RemoveSequences;
    for(auto &r: results) {
      report(cout, r, original_ds, options);
      res_t r2 = r;
      r2.counts = count_motif(ds, r.motif, options);
      report(cerr, r2, ds, options);
      cout << endl;
      Seeding::apply_mask(ds, r.motif, opts);
    }
  }

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

  return(EXIT_SUCCESS);
}

