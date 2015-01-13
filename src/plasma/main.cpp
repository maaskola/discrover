/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <sys/resource.h>  // for getrusage
#include <iostream>
#include <fstream>
#include <omp.h>
#include <boost/program_options.hpp>
#include "code.hpp"
#include "score.hpp"
#include "align.hpp"
#include "count.hpp"
#include "mask.hpp"
#include "plasma.hpp"
#include "../timer.hpp"
#include "cli.hpp"
#include "../executioninformation.hpp"
#include "../GitSHA1.hpp"
#include "../mcmc/montecarlo.hpp"

const std::string header = "# How to interpret this file:\n"
"# The program proceeds iteratively, at each step determining the single most discriminative word.";

const std::string program_name = "plasma";

std::string gen_usage_string() {
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

#if CAIRO_FOUND
Logo::matrix_t build_matrix(const string &motif, double absent) {
  const string nucls = "acgt";
  Logo::matrix_t matrix;
  for (auto pos : motif) {
    Logo::column_t col(4, 0);
    double z = 0;
    for (size_t i = 0; i < nucls.size(); i++)
      if (Seeding::iupac_included(nucls[i], pos)) {
        col[i] = 1;
        z++;
      }
    if (z == 0)
      for (size_t i = 0; i < nucls.size(); i++)
        col[i] = 1.0 / 4;
    else {
      for (size_t i = 0; i < nucls.size(); i++)
        if (col[i] > 0)
          col[i] = (1.0 - (4 - z) * absent) / z;
        else
          col[i] = absent;
    }
    matrix.push_back(col);
  }
  return matrix;
}

void generate_logos(const Seeding::Results &results,
                    const Seeding::Options &options) {
  size_t motif_idx = 0;
  for (auto &result : results) {
    Logo::matrix_t matrix = build_matrix(result.motif, options.logo.absent);
    Logo::draw_logo(matrix, options.label + ".motif" + to_string(motif_idx++),
                    options.logo);
  }
}
#endif

int main(int argc, const char **argv) {
  const string default_error_msg
      = "Please inspect the command line help with -h or --help.";
  Timer timer;

  Seeding::Options options;

  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Basic options");
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("verbose,v", "Be verbose about the progress")
    ("noisy,V", "Be very verbose about the progress")
    ;
  po::options_description ext_options = gen_plasma_options_description(options);

  desc.add(ext_options);

  po::positional_options_description pos;
  pos.add("fasta", -1);

  po::variables_map vm;

  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
  } catch (po::unknown_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " not known." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::ambiguous_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " is ambiguous." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_occurrences &e) {
    cout << "Error while parsing command line options:" << endl << "Option --"
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
         << "The value specified for option " << e.get_option_name()
         << " has an invalid format." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Too many positional options were specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
         << "Invalid command line syntax." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
         << "There is a programming error related to command line style."
         << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
         << "The configuration file can not be read." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Validation of option " << e.get_option_name() << " failed." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::error &e) {
    cout << "Error while parsing command line options:" << endl
         << "No further information as to the nature of this error is "
            "available, please check your command line arguments." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (exception &e) {
    cout << "An error occurred while parsing command line options." << endl
         << e.what() << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  if (vm.count("verbose"))
    options.verbosity = Verbosity::verbose;
  if (vm.count("noisy"))
    options.verbosity = Verbosity::debug;

  if (vm.count("version") and not vm.count("help")) {
    cout << program_name << " " << GIT_DESCRIPTION << " [" << GIT_BRANCH
         << " branch]" << endl;
    if (options.verbosity >= Verbosity::verbose)
      cout << GIT_SHA1 << endl;
    return EXIT_SUCCESS;
  }

  if (vm.count("help")) {
    cout << program_name << " " << GIT_DESCRIPTION << endl
         << "Copyright (C) 2011 Jonas Maaskola\n"
            "Provided under GNU General Public License Version 3 or later.\n"
            "See the file COPYING provided with this software for details of "
            "the license.\n" << endl;
    cout << gen_usage_string() << endl;
    cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  try {
    po::notify(vm);
  } catch (po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
         << "The value specified for option " << e.get_option_name()
         << " has an invalid format." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Too many positional options were specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
         << "Invalid command line syntax." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
         << "There is a programming error related to command line style."
         << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
         << "The configuration file can not be read." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::required_option &e) {
    cout << "Error while parsing command line options:" << endl
         << "The required option " << e.get_option_name()
         << " was not specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Validation of option " << e.get_option_name() << " failed." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::error &e) {
    cout << "Error while parsing command line options:" << endl
         << "No further information as to the nature of this error is "
            "available, please check your command line arguments." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  if (options.motif_specifications.size() == 0) {
    cout << "Error: you must specify motif lengths of interest with the -m "
            "switch." << endl;
    return EXIT_FAILURE;
  }

  if (vm.count("threads"))
    omp_set_num_threads(options.n_threads);

  if (not vm.count("output")) {
    options.label = generate_random_label(program_name, 0, options.verbosity);
    if (options.verbosity >= Verbosity::info)
      cout << "Using \"" << options.label
           << "\" as label to generate output file names." << endl;
  }

  if (options.verbosity >= Verbosity::verbose) {
    cout << "motif_specifications:";
    for (auto &x : options.motif_specifications)
      cout << " " << x;
    cout << endl;
    cout << "paths:";
    for (auto &x : options.paths)
      cout << " " << x;
    cout << endl;
    cout << "objectives:";
    for (auto &x : options.objectives)
      cout << " " << x;
    cout << endl;
  }

  try {
    // check and harmonize specified motifs, paths, and objectives
    Specification::harmonize(options.motif_specifications, options.paths,
                             options.objectives, false);
  } catch (runtime_error &e) {
    cout << e.what() << endl;
    return EXIT_FAILURE;
  }

  if (options.verbosity >= Verbosity::verbose) {
    cout << "motif_specifications:";
    for (auto &x : options.motif_specifications)
      cout << " " << x;
    cout << endl;
    cout << "paths:";
    for (auto &x : options.paths)
      cout << " " << x;
    cout << endl;
    cout << "objectives:";
    for (auto &x : options.objectives)
      cout << " " << x;
    cout << endl;
  }

  // initialize RNG
  if (options.verbosity >= Verbosity::verbose)
    cout << "Initializing random number generator with salt "
         << options.mcmc.random_salt << "." << endl;
  mt19937 rng;
  rng.seed(options.mcmc.random_salt);

  uniform_int_distribution<size_t> r_unif;

  Fasta::EntropySource::seed(r_unif(rng));
  MCMC::EntropySource::seed(r_unif(rng));

  Seeding::Plasma plasma(options);
  Seeding::Collection ds = plasma.collection;
  Seeding::Results results;
  using res_t = Seeding::Result;

  size_t n = options.motif_specifications.size();
  for (auto &motif : options.motif_specifications) {
    for (auto &r : plasma.find_motifs(
             motif, Seeding::objective_for_motif(options.objectives, motif),
             false))
      results.push_back(r);
    if (--n > 0)
      plasma.apply_mask(results);
  }

  if (results.size() == 1)
    report(cout, results[0], ds, options);
  else {
    sort(begin(results), end(results),
         [](const res_t &a, const res_t &b) { return a.log_p <= b.log_p; });
    Seeding::Collection original_ds = ds;
    Seeding::Options opts = options;
    opts.occurrence_filter = Seeding::OccurrenceFilter::RemoveSequences;
    for (auto &r : results) {
      report(cout, r, original_ds, options);
      res_t r2 = r;
      r2.counts = count_motif(ds, r.motif, options);
      report(cerr, r2, ds, options);
      cout << endl;
      try {
        Seeding::apply_mask(ds, r.motif, opts);
      } catch (exception &e) {
        cout << e.what() << endl;
        return EXIT_FAILURE;
      }
    }
  }

#if CAIRO_FOUND
  options.logo.revcomp = options.revcomp;
  generate_logos(results, options);
#endif

  if (options.verbosity >= Verbosity::info) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      cout << "getrusage failed" << endl;
      return EXIT_FAILURE;
    }

    double utime = usage.ru_utime.tv_sec + 1e-6 * usage.ru_utime.tv_usec;
    double stime = usage.ru_stime.tv_sec + 1e-6 * usage.ru_stime.tv_usec;
    double total_time = utime + stime;
    double elapsed_time = timer.tock() * 1e-6;

    cerr << "User time = " << utime << " sec" << endl
         << "System time = " << stime << " sec" << endl
         << "CPU time = " << total_time << " sec" << endl
         << "Elapsed time = " << elapsed_time << " sec" << endl
         << 100 * total_time / elapsed_time << "\% CPU" << endl;
  }
  return EXIT_SUCCESS;
}
