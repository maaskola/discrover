/* =====================================================================================
 * Copyright (c) 2011, Jonas Maaskola
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
 *       Filename:  main.cpp
 *
 *    Description:  Executable for the HMM package
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <sys/resource.h>  // for getrusage
#include <iostream>
#include <fstream>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../random_distributions.hpp"
#include "analysis.hpp"
#include "../aux.hpp"
#include "../terminal.hpp"
#include "../random_seed.hpp"
#include "../mcmc/montecarlo.hpp"
#include "../timer.hpp"
#include "../plasma/harmonization.hpp"
#include <git_config.hpp>
#include <discrover_paths.hpp>
#include "cli.hpp"

using namespace std;

string gen_usage_string(const string &program_name) {
  string usage = "This programs helps you discover binding site patterns in nucleic acids sequences with discriminative learning. "
    "Sets of positive and negative example sequences are mined for sequence motifs whose occurrence frequency varies between the sets. "
    "Hidden Markov models are used as probabilistic models of the sequence motifs. "
    "The method is applicable to the analysis of genome- and transcriptome-scale data, such as generated by next-generation sequencing technologies like ChIP-Seq or CLIP-Seq / PAR-CLIP. "
    "It makes use of available repeat experiments and aside from binary contrasts also more complex data configurations can be utilized. "
    "Several objective functions are available: for details on these please refer to the description of the --score option.\n"
    "\n"
    "The method is described in an open-access publication:\n"
    "Jonas Maaskola and Nikolaus Rajewsky. "
    "Binding site discovery from nucleic acid sequences by discriminative learning of hidden Markov models\n"
    "Nucleic Acid Research, 42(21):12995-13011, Dec 2014. doi:10.1093/nar/gku1083\n"
    "\n"
    "Please note that the software is accompanied by a manual in PDF format which includes recipes for motif discovery for ChIP-Seq and PAR-CLIP data.\n"
    "\n"
    "How to use this software.\n"
    "In the first example we want to train a motif seeded with the IUPAC string 'tgtanata' on sequences in the FASTA files signal.fa and control.fa. "
    // "By default, it will train in hybrid mode, using mutual information for the emission probabilities of the motifs, and likelihood maximization for everything else (background emission and all transition probabilities). "
    "Output will be written to files whose file name is prefixed by 'label':\n"
    "\n"
    "" + program_name + " signal.fa control.fa -m tgtanata -o label\n"
    "\n"
    "Similarly, to automatically find a motif of length 8 one can use the following:\n"
    "\n"
    "" + program_name + " signal.fa control.fa -m 8 -o label\n"
    "\n"
    "When the option --output is not given, a unique output prefix will be automatically generated.\n"
    "\n"
    "It is possible to use just a single FASTA file for discriminative sequence analysis. "
    "In this case a control set of sequences will be generated by shuffling the signal sequences:\n"
    "\n"
    "" + program_name + " signal.fa -m 8\n"
    "\n"
     "Note that names may be given to the motifs to annotate in which of the samples they are "
    "respectively expected to be enriched, as the following example demonstrates.\n"
    "\n"
    "" + program_name + " A:sample1.fa B:sample2.fa -m A:8 -m B:6\n"
    "\n"
    "Here, A and B are arbitrary labels given to motifs of length 8 and 6 that are to be enriched in sample1 and sample2, respectively.\n"
    "\n"
    "Options may also be specified in a configuration file, by using --config. "
    "Options given as arguments on the command line take precedence over configuration file entries. "
    "For options accepted multiple times, parameters from command line and configuration file are joined.\n"
    "\n"
    "Various output files are generated, including the resulting parameters, summary information, motif occurrence tables, and the Viterbi parse. "
#if CAIRO_FOUND
    "Also, sequence logos of the found motifs are automatically generated. "
#endif
    "See the description of the --output option below for details.\n"
    "\n"
    "Apart from the description of command line arguments below, advanced usage instructions are given in the manual. "
    "In particular, the manual explains\n"
    "a) how sets of sequences can be grouped into multiple contrasts (useful e.g. when repeat experiments are available),\n"
    "b) how motifs can be discovered that are discriminative for one contrast but not for another, and\n"
    "c) how composite models of multiple motifs can be built, which is useful to discover co-factor motifs e.g. in ChIP-Seq analysis (in brief: by discovering many seeds and then building composite HMMs using the --multiple switch).\n"
    "\n"
    "Manual: "
    DISCROVER_MANUAL_PATH
    ;
  return usage;
}

Measures::Discrete::Measure measure2iupac_objective(const Measure measure) {
  if (Measures::is_generative(measure))
    return Measures::Discrete::Measure::SignalFrequency;
  switch (measure) {
    case Measure::DeltaFrequency:
      return Measures::Discrete::Measure::DeltaFrequency;
    case Measure::MatthewsCorrelationCoefficient:
      return Measures::Discrete::Measure::MatthewsCorrelationCoefficient;
    default:
      return Measures::Discrete::Measure::MutualInformation;
  }
}

void fixup_seeding_options(Options::HMM &options) {
  // TODO REACTIVATE if(set_objective)
  // TODO REACTIVATE    options.seeding.objective =
  // {measure2iupac_objective(options.measure);
  // options.seeding.objective = Seeding::Objective::corrected_logp_gtest;
  if (options.seeding.objectives.size() == 1
      and begin(options.seeding.objectives)->measure
          == Measures::Discrete::Measure::SignalFrequency)
    options.seeding.plasma.rel_degeneracy = 0.2;
  options.seeding.paths = options.paths;
  options.seeding.n_threads = options.n_threads;
  options.seeding.n_seq = options.n_seq;
  options.seeding.weighting = options.weighting;
  // options.seeding.verbosity = Verbosity::info;
  options.seeding.verbosity = options.verbosity;
  options.seeding.revcomp = options.revcomp;
  options.seeding.pseudo_count = options.contingency_pseudo_count;
  options.seeding.measure_runtime = options.timing_information;
  options.seeding.label = options.label;
#if CAIRO_FOUND
  options.logo.revcomp = options.revcomp;
  options.seeding.logo = options.logo;
  options.seeding.logo.pdf_logo = false;
  options.seeding.logo.png_logo = false;
#endif
  options.seeding.mcmc.max_iter = options.termination.max_iter;
  options.seeding.mcmc.temperature = options.sampling.temperature;
  options.seeding.mcmc.n_parallel = options.sampling.n_parallel;
  options.seeding.mcmc.random_salt = options.random_salt;
}

int main(int argc, const char **argv) {
  const string default_error_msg
      = "Please inspect the command line help with -h or --help.";
  Timer timer;

  namespace po = boost::program_options;

  Options::HMM options;
  options.n_threads = omp_get_num_procs();
  options.exec_info
      = ExecutionInformation(argv[0], GIT_DESCRIPTION, GIT_BRANCH, argc, argv);
  options.class_model = false;
  options.random_salt = generate_rng_seed();

  string config_path;

  string background_initialization;

  static const size_t MIN_COLS = 60;
  static const size_t MAX_COLS = 80;
  size_t cols = get_terminal_width();
  if (cols < MIN_COLS)
    cols = MIN_COLS;
  if (cols > MAX_COLS)
    cols = MAX_COLS;

  po::options_description hidden_options("Hidden options", cols);
  po::options_description common_options, visible_options, cmdline_options,
      config_file_options;
  try {
    gen_discrover_cli(cols, config_path, options, common_options,
                      visible_options, hidden_options, cmdline_options,
                      config_file_options);
  } catch (...) {
    cout << "Error while generating command line options." << endl
         << "Please notify the developers." << endl;
    return EXIT_FAILURE;
  }

  po::positional_options_description pos;
  pos.add("fasta", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(cmdline_options)
                  .positional(pos)
                  .run(),
              vm);
  } catch (po::unknown_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " not known." << endl << default_error_msg
         << endl;
    return EXIT_FAILURE;
  } catch (po::ambiguous_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " is ambiguous." << endl << default_error_msg
         << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_occurrences &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
         << "The value specified for option " << e.get_option_name()
         << " has an invalid format." << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Too many positional options were specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
         << "Invalid command line syntax." << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
         << "There is a programming error related to command line style."
         << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
         << "The config file can not be read." << endl << default_error_msg
         << endl;
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

  options.verbosity = Verbosity::info;
  if (vm.count("verbose"))
    options.verbosity = Verbosity::verbose;
  if (vm.count("noisy"))
    options.verbosity = Verbosity::debug;

  if (vm.count("version") and not vm.count("help")) {
    cout << options.exec_info.program_name << " "
         << options.exec_info.hmm_version << " [" << GIT_BRANCH << " branch]"
         << endl;
    if (options.verbosity >= Verbosity::verbose)
      cout << GIT_SHA1 << endl;
    return EXIT_SUCCESS;
  }

  if (vm.count("help")) {
    cout << options.exec_info.program_name << " "
         << options.exec_info.hmm_version << endl;
    cout << "Copyright (C) 2011 Jonas Maaskola\n"
            "Provided under GNU General Public License Version 3 or later.\n"
            "See the file COPYING provided with this software for details of "
            "the license.\n" << endl;
    cout << limit_line_length(gen_usage_string(options.exec_info.program_name),
                              cols) << endl << endl;
    cout << visible_options << endl;
    switch (options.verbosity) {
      case Verbosity::nothing:
      case Verbosity::error:
      case Verbosity::info:
        cout << "Advanced and hidden options not shown. Use -hv or -hV to show "
                "them." << endl;
        break;
      case Verbosity::verbose:
        cout << common_options << endl;
        cout << "Hidden options not shown. Use -hV to show them." << endl;
        break;
      case Verbosity::debug:
      case Verbosity::everything:
        cout << common_options << endl;
        cout << hidden_options << endl;
        break;
    }
    return EXIT_SUCCESS;
  }

  try {
    po::notify(vm);
  } catch (po::required_option &e) {
    cout << "Error while parsing command line options:" << endl
         << "The required option " << e.get_option_name()
         << " was not specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  if (config_path != "") {
    ifstream ifs(config_path.c_str());
    if (!ifs) {
      cout << "Error: can not open config file: " << config_path << endl
           << default_error_msg << endl;
      return EXIT_FAILURE;
    } else {
      try {
        store(parse_config_file(ifs, config_file_options), vm);
      } catch (po::multiple_occurrences &e) {
        cout << "Error while parsing config file:" << endl << "Option "
             << e.get_option_name() << " was specified multiple times." << endl
             << default_error_msg << endl;
        return EXIT_FAILURE;
      } catch (po::unknown_option &e) {
        cout << "Error while parsing config file:" << endl << "Option "
             << e.get_option_name() << " not known." << endl
             << default_error_msg << endl;
        return EXIT_FAILURE;
      } catch (po::invalid_option_value &e) {
        cout << "Error while parsing config file:" << endl
             << "The value specified for option " << e.get_option_name()
             << " has an invalid format." << endl
             << default_error_msg << endl;
        return EXIT_FAILURE;
      }
      notify(vm);
    }
  }

  // generate an output path stem if the user did not specify one
  if (not vm.count("output")) {
    options.label = generate_random_label(options.exec_info.program_name, 0,
                                          options.verbosity);
    while (boost::filesystem::exists(options.label + ".hmm"))
      options.label = generate_random_label(options.exec_info.program_name, 5,
                                            options.verbosity);
    if (options.verbosity >= Verbosity::info)
      cout << "Using \"" << options.label
           << "\" as label to generate output file names." << endl;
  }

  // set the number of threads with OpenMP
  omp_set_num_threads(options.n_threads);

  // print information about specified motifs, paths, and objectives
  if (options.verbosity >= Verbosity::debug) {
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
                             options.objectives);
  } catch (runtime_error &e) {
    cout << e.what() << endl;
    return EXIT_FAILURE;
  }

  // print information about specified motifs, paths, and objectives
  if (options.verbosity >= Verbosity::debug) {
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

  // ensure a positive number, if anything, is specified for sampling min_size
  if (vm.count("smin")) {
    if (options.sampling.min_size < 0) {
      cout << "Please note that the minimal length for MCMC sampling must be "
              "non-negative." << endl << default_error_msg << endl;
      return EXIT_FAILURE;
    }
  } else
    options.sampling.min_size = -1;

  // ensure a positive number, if anything, is specified for sampling max_size
  if (vm.count("smax")) {
    if (options.sampling.max_size < 0) {
      cout << "Please note that the maximal length for MCMC sampling must be "
              "non-negative." << endl << default_error_msg << endl;
      return EXIT_FAILURE;
    } else if (options.sampling.max_size < options.sampling.min_size) {
      cout << "Please note that the maximal length for MCMC sampling must be "
              "larger than the minimal size." << endl
              << default_error_msg << endl;
      return EXIT_FAILURE;
    }
  } else
    options.sampling.max_size = -1;

  // Initialize the plasma options
  fixup_seeding_options(options);

  if (options.termination.past == 0) {
    cout << "Error: the value of --past must be a number greater than 0."
         << endl;
  }

  if (options.termination.max_iter == 0 and options.sampling.do_sampling) {
    options.termination.max_iter = 1000;
    cout << "Note: did not specify the number of iterations to perform "
            "(--maxiter)." << endl << "We will now do "
         << options.termination.max_iter << " iterations." << endl;
  }

  bool any_named = false;
  for (auto &x : options.paths)
    if (not x.motifs.empty()) {
      any_named = true;
      break;
    }
  if (not any_named)
    for (auto &x : options.paths)
      for (auto &s : options.motif_specifications)
        x.motifs.insert(s.name);

  if (options.load_paths.empty() and options.motif_specifications.empty()) {
    cout << "Error: you must either specify at least one of:" << endl
         << "1. a path from which to load HMM parameter (--load)" << endl
         << "2. one or more motifs (--motif)" << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  // Ensure that the residual MI ratio cutoff is non-negative
  if (options.multi_motif.residual_ratio < 0) {
    cout << "Warning: negative value provided for residual mutual information "
            "ratio cutoff. Using 0 as value." << endl;
    options.multi_motif.residual_ratio = 0;
  }

  // Ensure that multiple mode is only used with objective function MICO
  if (options.multi_motif.accept_multiple) {
    for (auto &obj : options.objectives)
      if (obj.measure != Measures::Continuous::Measure::MutualInformation) {
        cout << "Error: multiple motif mode can only be used with the "
                "objective function MICO." << endl;
        return EXIT_FAILURE;
      }
  }

  if (options.line_search.eta <= options.line_search.mu) {
    cout << "Error: the Moré-Thuente η parameter must be larger than the µ "
            "parameter." << endl;
    return EXIT_FAILURE;
  }

  // initialize RNG
  if (options.verbosity >= Verbosity::info)
    cout << "Initializing random number generator with salt "
         << options.random_salt << "." << endl;
  mt19937 rng;
  rng.seed(options.random_salt);

  Fasta::EntropySource::seed(RandomDistribution::Uniform(rng));
  MCMC::EntropySource::seed(RandomDistribution::Uniform(rng));

  // main routine
  try {
    perform_analysis(options, rng);
  } catch (exception &e) {
    cout << e.what() << endl;
    return EXIT_FAILURE;
  }

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

    const size_t label_col_width = 15;
    const size_t value_col_width = 15;
    cerr << setw(label_col_width) << left << "User time"
         << setw(value_col_width) << right << utime << " sec"
         << endl
         << setw(label_col_width) << left << "System time"
         << setw(value_col_width) << right << stime << " sec"
         << endl
         << setw(label_col_width) << left << "CPU time"
         << setw(value_col_width) << right << total_time << " sec"
         << endl
         << setw(label_col_width) << left << "Elapsed time"
         << setw(value_col_width) << right << elapsed_time << " sec"
         << endl
         << setw(label_col_width) << left << "CPU \%"
         << setw(value_col_width) << right << 100 * total_time / elapsed_time
         << endl;
  }

  return EXIT_SUCCESS;
}
