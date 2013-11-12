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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <ctime>          // for clock()
#include <sys/times.h>    // for times(struct tms *)
#include <iostream>
#include <fstream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "plasma/plasma_cli.hpp"
// #include "plasma/GitSHA1.hpp"
#include "GitSHA1.hpp"
#include "analysis.hpp"
#include "aux.hpp"
#include "terminal.hpp"

using namespace std;

string gen_usage_string(const string &program_name)
{
  string usage = "This program implements hidden Markov models for probabilistic sequence analysis, "
    "in particular for the purpose of discovering unknown binding site patterns in nucleic acid sequences. "
    "It may be used to train models using sequence data, evaluate models on "
    "sequence data, or to simulate sequence data from a model. "
    "Several learning objectives are implemented. "
    "Please refer to the description for the --score option below.\n"
   "\n"
    "Among these objectives there are two non-discriminative measures: "
    "Viterbi learning and maximum likelihood learning using the Baum-Welch method. "
    "The other measures are discriminative and use gradient-based learning.\n"
    "The discriminative measures require the specification of multiple data sets. "
    "\n"
//    "Except for mutual information, the discriminative methods require paired signal and control "
//    "data sets.\n"
    "Here is an example call to train a motif seeded with the IUPAC string 'tgtanata' on "
    "sequences in files signal.fa and control.fa. By default, it will train in hybrid mode, using "
    "mutual information for the emission probabilities of the motifs, and likelihood maximization "
    "for everything else (background emission and all transition probabilities). Output will be "
    "written to files whose file name is prefixed by 'output':\n"
    "\n"
    "" + program_name + " -f signal.fa -f control.fa -m tgtanata -o output\n"
    "\n"
    "Similarly, to automatically find a motif of length 8 one can use the following:\n"
    "\n"
    "" + program_name + " -f signal.fa -f control.fa -m 8 -o output\n"
    "\n"
    "Actually, the -f switches are not needed, as any free arguments will be treated as paths to FASTA files. Consequently, the same command may be given as:\n"
    "\n"
    "" + program_name + " signal.fa control.fa -m 8 -o output\n"
    "\n"
    "Note that you may give names to the motifs, and annotate in which of the samples they are "
    "respectively expected to be enriched, as the following example demonstrates.\n"
    "\n"
    "" + program_name + " pum_motif:pum.signal.fa pum.control.fa qki_motif:qki.signal.fa qki.control.fa -m pum_motif:tgtanata -m qki_motif:nnntaacn -o output\n"
    "\n"
    "Here, pum_motif and qki_motif are arbitrary labels given to the motifs. In the example two of the FASTA sequence files are expected to be enriched for the respective motif. "
    "Note that it is necessary to annotate sequence files in this manner in order to evaluate all discriminative measures except for mutual information, or to train using these as objective functions.\n"
    "\n"
    "Options may also be specfied in a configuration file, by using --config. "
    "Options given as arguments on command line take precedence over configuration "
    "file entries. For options accepted multiple times parameters from command line "
    "and configuration file are joined.\n"
    "\n"
    "By default, all available CPU cores will be used by starting as many threads as there are "
    "CPU cores on the machine. To start a specific number of threads the -t switch or the "
    "environment variable OMP_NUM_THREADS may be used. The following examples use four threads:\n"
    "\n"
    "" + program_name + " signal.fa control.fa -m tgtanata -o output -t 4\n"
    "OMP_NUM_THREADS=4 " + program_name + " signal.fa control.fa -m tgtanata -o output\n"
    "\n"
    "Various output files are generated, including the resulting parameters, summary information, motif occurrence tables, and the Viterbi parse. See the description of the required argument -o / --output below for details.";
  return(usage);
}

Measures::Discrete::Measure measure2iupac_objective(const Measure measure) {
  if(Measures::is_generative(measure))
    return(Measures::Discrete::Measure::SignalFrequency);
  switch(measure) {
    case Measure::DeltaFrequency:
      return(Measures::Discrete::Measure::DeltaFrequency);
    case Measure::MatthewsCorrelationCoefficient:
      return(Measures::Discrete::Measure::MatthewsCorrelationCoefficient);
    default:
      return(Measures::Discrete::Measure::MutualInformation);
  }
}

void fixup_seeding_options(hmm_options &options) {
// TODO REACTIVATE if(set_objective)
// TODO REACTIVATE    options.seeding.objective = {measure2iupac_objective(options.measure);
  // options.seeding.objective = Seeding::Objective::corrected_logp_gtest;
  if(options.seeding.objectives.size() == 1 and begin(options.seeding.objectives)->measure == Measures::Discrete::Measure::SignalFrequency)
    options.seeding.rel_degeneracy = 0.2;
  options.seeding.paths = options.paths;
  options.seeding.n_threads = options.n_threads;
  options.seeding.n_seq = options.n_seq;
  // options.seeding.verbosity = Verbosity::info;
  options.seeding.verbosity = options.verbosity;
  options.seeding.revcomp = options.revcomp;
  options.seeding.pseudo_count = options.contingency_pseudo_count;
  options.seeding.measure_runtime = options.timing_information;
}

string generate_random_label(const string &prefix="dlhmm", size_t n_rnd_char=5, Verbosity verbosity=Verbosity::info) {
  using namespace boost::posix_time;
  ptime t = microsec_clock::universal_time();
  string datetime = to_iso_extended_string(t) + "Z";

  // TODO perhaps add the process ID -> getpid()
  if(verbosity >= Verbosity::debug)
    cout << "Generating random label with prefix " << prefix << " and " << n_rnd_char << " random characters." << endl;
  std::string label = prefix + "_" + datetime;
  if(n_rnd_char > 0) {
    label += "_";
    const char first = 'a';
    const char last = 'z';
    size_t n = last - first + 1;
    for(size_t i = 0; i < n_rnd_char; i++)
      label += char(rand() % n + first);
  }
  if(verbosity >= Verbosity::debug)
    cout << "Generated random label " << label << endl;
  return(label);
}

template <typename X> X entropy_from_os(const std::string &source="/dev/urandom") {
  std::ifstream f(source.c_str());
  X myRandom;
  f.read(reinterpret_cast<char*>(&myRandom), sizeof(myRandom));
  return(myRandom);
}

template <typename X> void mix_with_os_entropy(X &x, const std::string &source="/dev/urandom") {
  if(boost::filesystem::exists(source)) {
    X y = entropy_from_os<X>(source);
    x = x ^ y;
  }
}

int main(int argc, const char** argv)
{
  const string default_error_msg = "Please inspect the command line help with -h or --help.";
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  namespace po = boost::program_options;

  hmm_options options;
  options.exec_info = generate_exec_info(argv[0], GIT_DESCRIPTION, cmdline(argc, argv));
  options.class_model = false;
  options.random_salt = time(0) ^ getpid(); // XOR Unix time & process ID
  mix_with_os_entropy(options.random_salt); // XOR with entropy from /dev/urandom if available

  string config_path;
  bool hmm_score_seeding;

  string background_initialization;

  static const size_t MIN_COLS = 60;
  static const size_t MAX_COLS = 100;
  size_t cols = get_terminal_width();
  if(cols < MIN_COLS)
    cols = MIN_COLS;
  if(cols > MAX_COLS)
    cols = MAX_COLS;

   // Declare the supported options.
  po::options_description generic_options("Generic options", cols);
  po::options_description basic_options("Basic options, required", cols);
  po::options_description objective_options("Objective function options", cols);
  po::options_description init_options("Initialization options", cols);
  po::options_description eval_options("Evaluation options", cols);
  po::options_description advanced_options("Advanced options", cols);
  po::options_description sampling_options("Gibbs sampling options", cols);
  po::options_description termination_options("Termination options", cols);
  po::options_description hidden_options("Hidden options", cols);

  po::options_description seeding_options = gen_iupac_options_description(options.seeding, "", "Seeding options for IUPAC regular expression finding", cols, false, false);

  generic_options.add_options()
    ("config", po::value<string>(&config_path), "Read options from a configuration file. ")
    ("help,h", "Produce help message. Combine with -v or -V for additional commands.")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("verbose,v", "Be verbose about the progress.")
    ("noisy,V", "Be very verbose about the progress.")
    ;

  basic_options.add_options()
    ("fasta,f", po::value<vector<Specification::DataSet>>(&options.paths)->required(), "Path to FASTA file with sequences. May be given multiple times. The path may be prepended by a comma separated list of names of motifs enriched in the file, with the syntax [NAMES:[SERIES:]]PATH, where NAMES is NAME[,NAMES] is a set of motif names, and SERIES is the name of a series this data set belongs to. Also see the example above. Note that you should prepend the path with a double colon if the path contains at least one colon.") // NOTE: one might not require this argument in case only simulation is wanted, but for the usual case this is required
    ("output,o", po::value<string>(&options.label), "Output file names are generated from this label. If this option is not specified the output label will be 'dlhmm_XXX' where XXX is a random string to make the label unique. The output files comprise:\n"
     ".hmm     \tParameter of the trained HMM. May be loaded later with -l.\n"
     ".summary \tSummary information with number of occurrences of the motifs in each sample, and various generative and discriminative statistics.\n"
     ".viterbi \tFASTA sequences annotated with the viterbi path, and sequence level statistics.\n"
     ".table   \tCoordinates and sequences of matches to the motifs in all sequences.\n"
     "Note that, depending on the argument of --compress, the latter two files may be compressed, and require decompression for inspection.")
    ("motif,m", po::value<vector<Specification::Motif>>(&options.motif_specifications), "Motif specification. May be given multiple times, and can be specified in multiple ways:\n"
        "1. \tUsing the IUPAC code for nucleic acids.\n"
        "2. \tA length specification. Motifs of the indicated lengths are sought. A length specification is a comma separated list of length ranges, where a length range is either a single length, or an expression of the form 'x-y' to indicate lengths x up to y.\n"
        "3. \tBy specifying the path to a file with a PWM.\n"
        "Regardless of the way a motif is specified, it may be given a name, and an insertions string. The syntax is [NAME:[INSERT:]]MOTIFSPEC. If MOTIFSPEC is a path and contains at least one colon please just prepend two colons. The insertions string is a comma separated list of positions after which to add an insertion state. Positions are 1 indexed, and must be greater 1 and less than the motif length.")
    ;

  objective_options.add_options()
    ("score", po::value<Training::Objectives>(&options.objectives)->default_value(Training::Objectives(1,Training::Objective("mi")), "mi"), "The significance measure. May be one of\n"
//  "likelihood \tLikelihood using Baum-Welch\n"
//   "cml     \tCorrect classification probability\n"
//   "chisq   \tChi-square\n"
     "none    \tDo not perform learning\n"
     "bw      \tLikelihood using Baum-Welch\n"
     "viterbi \tViterbi learning\n"
     "mi      \tMutual information of condition and motif occurrence\n"
     "ri      \tRank mutual information\n"
     "mmie    \tMaximum mutual information (MMIE): posterior classification probability\n"
     "mcc     \tMatthew's correlation coefficient\n"
     "dlogl   \tLog-likelihood difference, like DME, see doi: 10.1073/pnas.0406123102\n"
     "dfreq   \tDifference of frequency of sequences with motif occurrences, similar to DIPS and DECOD, see doi: 10.1093/bioinformatics/btl227 and 10.1093/bioinformatics/btr412\n")
//    ("class", po::bool_switch(&options.class_model), "Whether to restrict to the annotated motifs during Baum-Welch updates.")
    ("revcomp,r", po::bool_switch(&options.revcomp), "Respect also the reverse complementary sequence.")
    ;

  init_options.add_options()
    ("load,l", po::value<vector<string>>(&options.load_paths), "Load HMM parameters from a .hmm file produced by an earlier run. Can be specified multiple times; then the first parameter file will be loaded, and motifs of the following parameter files are added.")
    ("alpha", po::value<double>(&options.alpha)->default_value(0.03, "0.03"), "Probability of alternative nucleotides. The nucleotides not included in the IUPAC character will have this probability.")
    ("lambda", po::value<double>(&options.lambda)->default_value(1), "Initial value for prior with which a motif is expected.")
    ("bgorder", po::value<size_t>(&options.bg_order)->default_value(0), "Order of the background model.")
    ("wiggle", po::value<size_t>(&options.wiggle)->default_value(0), "For automatically determined seeds, consider variants shifted up and down by up to the specified number of positions.")
    ("padl", po::value<size_t>(&options.left_padding)->default_value(0), "For automatically determined seeds, add Ns upstream of, or to the left of the seed.")
    ("padr", po::value<size_t>(&options.right_padding)->default_value(0), "For automatically determined seeds, add Ns downstream of, or to the right of the seed.")
    ("hmmscore", po::bool_switch(&hmm_score_seeding), "Train a HMM for each seed and select seeds in order of their HMM scores.")
    ("miseeding", po::bool_switch(&options.use_mi_to_seed), "Disregard automatic seeding choice and use MICO for seeding.")
    ;

  eval_options.add_options()
    ("no-summary", po::bool_switch(&options.evaluate.summary)->default_value(true), "Do not print summary information.")
    ("no-viterbi", po::bool_switch(&options.evaluate.viterbi_path)->default_value(true), "Do not print the Viterbi path.")
    ("no-occurrence-table", po::bool_switch(&options.evaluate.occurrence_table)->default_value(true), "Do not print the occurrence table.")
    ("no-ric", po::bool_switch(&options.evaluate.ric)->default_value(true), "Do not print the RIC analysis.")
    ;


  advanced_options.add_options()
    ("bg_learn", po::value<Training::Method>(&options.bg_learning)->default_value(Training::Method::Reestimation, "em"), "How to learn the background. Available are 'fixed', 'em', 'gradient', where the 'em' uses reestimation to maximize the likelihood contribution of the background parameters, while 'gradient' uses the discriminative objective functon.")
    ("ls_mu", po::value<double>(&options.line_search.mu)->default_value(0.1, "0.1"), "The parameter µ for the Moré-Thuente line search algorithm.")
    ("ls_eta", po::value<double>(&options.line_search.eta)->default_value(0.5, "0.5"), "The parameter η for the Moré-Thuente line search algorithm.")
    ("ls_delta", po::value<double>(&options.line_search.delta)->default_value(0.66, "0.66"), "The parameter delta for the Moré-Thuente line search algorithm.")
    ("ls_num", po::value<size_t>(&options.line_search.max_steps)->default_value(10, "10"), "How many gradient and function evaluation to perform maximally per line search.")
    ("mic", po::value<size_t>(&options.mic)->default_value(0), "How many additional rank groups to introduce in the data. Allows to optimize the maximal information content (MIC).")
    ("multi", po::value<Training::Simultaneity>(&options.simultaneity)->default_value(Training::Simultaneity::Simultaneous,"sim"), "Wether to add and train automatically determined motifs sequentially or simultaneously. Available are 'seq', 'sim'.")
    ("threads,T", po::value<size_t>(&options.n_threads)->default_value(omp_get_num_procs()), "The number of threads to use. If this in not specified, the value of the environment variable OMP_NUM_THREADS is used if that is defined, otherwise it will use as many as there are CPU cores on this machine.")
    ("runtime", po::bool_switch(&options.timing_information), "Output information about how long certain parts take to execute.")
    ("simulate", po::value<size_t>(&options.n_simulations)->default_value(0), "Simulate the model.")
    ("posterior", po::bool_switch(&options.print_posterior), "During evaluation also print out the motif posterior probabilitity.")
    ("no-classp", po::bool_switch(&options.learn_class_prior)->default_value(true), "When performing MMIE, do not learn the class prior.")
    ("no-motifp", po::bool_switch(&options.learn_conditional_motif_prior)->default_value(true), "When performing MMIE, do not learn the conditional motif prior.")
    ("classp", po::value<double>(&options.class_prior)->default_value(0.5), "When performing MMIE, use this as initial class prior.")
    ("motifp1", po::value<double>(&options.conditional_motif_prior1)->default_value(0.6, "0.6"), "When performing MMIE, use this as initial conditional motif prior for the signal class.")
    ("motifp2", po::value<double>(&options.conditional_motif_prior2)->default_value(0.03, "0.03"), "When performing MMIE, use this as initial conditional motif prior for the control class.")
    ("cv", po::value<size_t>(&options.cross_validation_iterations)->default_value(0), "Number of cross validation iterations to do.")
    ("cv_freq", po::value<double>(&options.cross_validation_freq)->default_value(0.9, "0.9"), "Fraction of data samples to use for training in cross validation.")
    ("nseq", po::value<size_t>(&options.n_seq)->default_value(0), "Only consider the top sequences of each FASTA file. Specify 0 to indicate all.")
    ("compress", po::value<Compression>(&options.output_compression)->default_value(Compression::gzip, "gz"), "Compression method to use for larger output files. Available are: 'none', 'gz' or 'gzip', 'bz2' or 'bzip2'.") // TODO make the code conditional on the presence of zlib
    ("pscnt", po::value<double>(&options.contingency_pseudo_count)->default_value(1.0, "1"), "The pseudo count to be added to the contingency tables in the discriminative algorithms.")
    ("pscntE", po::value<double>(&options.emission_pseudo_count)->default_value(1.0, "1"), "The pseudo count to be added to the expected emission probabilities before normalization in the Baum-Welch algorithm.")
    ("pscntT", po::value<double>(&options.transition_pseudo_count)->default_value(0.0, "0"), "The pseudo count to be added to the expected transition probabilities before normalization in the Baum-Welch algorithm.")
    ;

  termination_options.add_options()
    ("maxiter", po::value<size_t>(&options.termination.max_iter)->default_value(1000), "Maximal number of iterations to perform in training. A value of 0 means no limit, and that the training is only terminated by the tolerance.")
    ("gamma", po::value<double>(&options.termination.gamma_tolerance)->default_value(1e-4, "1e-4"), "Tolerance for the reestimation type learning methods. Training stops when the L1 norm of the parameter change between iterations is less than this value.")
    ("delta", po::value<double>(&options.termination.delta_tolerance)->default_value(1e-4, "1e-4"), "Relative score difference criterion tolerance for training algorithm termination: stops iterations when (f - f') / f < delta, where f' is the objective value of the past iteration, and f is the objective value of the current iteration.")
    ("epsilon", po::value<double>(&options.termination.epsilon_tolerance)->default_value(0), "Gradient norm criterion tolerance for training algorithm termination: stops when ||g|| < epsilon * max(1, ||g||), where ||.|| denotes the Euclidean (L2) norm.")
    ("past", po::value<size_t>(&options.termination.past)->default_value(1), "Distance for delta-based convergence test. This parameter determines the distance, in iterations, to compute the rate of decrease of the objective function.")
    ;

  sampling_options.add_options()
    ("sampling", po::bool_switch(&options.sampling.do_sampling), "Perform Gibbs sampling for parameter inference instead of re-estimation or gradient learning.")
    ("temp", po::value<double>(&options.sampling.temperature)->default_value(1e-3), "When performing Gibbs sampling use this temperature. The temperatures of parallel chains is decreasing by factors of two.")
//    ("anneal", po::value<double>(&options.sampling.anneal_factor)->default_value(0.98,"0.98"), "When performing Gibbs sampling multiply the temperature by this term in each iteration.")
    ("smin", po::value<int>(&options.sampling.min_size), "Minimal length for Gibbs sampling. When unspecified defaults to initial motif length.")
    ("smax", po::value<int>(&options.sampling.max_size), "Maximal length for Gibbs sampling. When unspecified defaults to initial motif length.")
    ("nindel", po::value<size_t>(&options.sampling.n_indels)->default_value(5), "Maximal number of positions that may be added or removed at a time. Adding and removing of happens at and from the ends of the motif.")
    ("nshift", po::value<size_t>(&options.sampling.n_shift)->default_value(5), "Maximal number of positions that the motif may be shifted by.")
    ("partemp", po::value<size_t>(&options.sampling.n_parallel)->default_value(6), "Parallel chains to run for parallel tempering.")
    ;

  hidden_options.add_options()
    ("absthresh", po::bool_switch(&options.termination.absolute_improvement), "Whether improvement should be gauged by absolute value. Default is relative to the current score.")
    ("salt", po::value<unsigned int>(&options.random_salt), "Seed for the random number generator. If unspecified, will use current time, XORed with the process ID and additionally with entropy from /dev/urandom, if available.")
    ("intermediate", po::bool_switch(&options.store_intermediate), "Write out intermediate parameters during training.")
    ("limitlogp", po::bool_switch(&options.limit_logp), "Do not report corrected log-P values greater 0 but report 0 in this case.")
    ("longnames", po::bool_switch(&options.long_names), "Form longer output file names that contain some information about the parameters.")
    ;

  advanced_options.add(termination_options);
  advanced_options.add(sampling_options);

  po::options_description simple_options;
  simple_options.add(basic_options).add(objective_options).add(init_options).add(seeding_options).add(eval_options);

  po::options_description common_options;
  common_options.add(simple_options).add(advanced_options);

  po::options_description visible_options;
  visible_options.add(generic_options).add(simple_options);

  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(common_options).add(hidden_options);

  po::options_description config_file_options;
  config_file_options.add(common_options).add(hidden_options);


  po::positional_options_description pos;
  pos.add("fasta", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  } catch(po::unknown_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " not known." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::ambiguous_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " is ambiguous." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " was specified multiple times." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::multiple_occurrences &e) {
    cout << "Error while parsing command line options:" << endl
      << "Option " << e.get_option_name() << " was specified multiple times." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
      << "The value specified for option " << e.get_option_name() << " has an invalid format." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Too many positional options were specified." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
      << "Invalid command line syntax." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
      << "There is a programming error related to command line style." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
      << "The config file can not be read." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
      << "Validation of option " << e.get_option_name() << " failed." << endl
      << default_error_msg << endl;
    return(-1);
  } catch(po::error &e) {
    cout << "Error while parsing command line options:" << endl
      << "No further information as to the nature of this error is available, please check your command line arguments." << endl
      << default_error_msg << endl;
    return(-1);
  }

  options.verbosity = Verbosity::info;
  if(vm.count("verbose"))
    options.verbosity = Verbosity::verbose;
  if(vm.count("noisy"))
    options.verbosity = Verbosity::debug;

  if(vm.count("version") and not vm.count("help")) {
    cout << options.exec_info.program_name << " " << options.exec_info.hmm_version << endl;
    if(options.verbosity >= Verbosity::verbose)
      cout << GIT_SHA1 << endl;
    return(EXIT_SUCCESS);
  }

  if(vm.count("help")) {
    cout << options.exec_info.program_name << " " << options.exec_info.hmm_version << endl;
    cout << "Copyright (C) 2011 Jonas Maaskola\n"
      "Provided under GNU General Public License Version 3 or later.\n" 
      "See the file COPYING provided with this software for details of the license.\n" << endl;
    cout << limit_line_length(gen_usage_string(options.exec_info.program_name), cols) << endl << endl;
    cout << visible_options << endl;
    switch(options.verbosity) {
      case Verbosity::nothing:
      case Verbosity::error:
      case Verbosity::info:
        cout << "Advanced and hidden options not shown. Use -hv or -hV to show them." << endl;
        break;
      case Verbosity::verbose:
        cout << advanced_options << endl;
        cout << "Hidden options not shown. Use -hV to show them." << endl;
        break;
      case Verbosity::debug:
      case Verbosity::everything:
        cout << advanced_options << endl;
        cout << hidden_options << endl;
        break;
    }
    return 1;
  }

  try {
    po::notify(vm);
  } catch(po::required_option &e) {
    cout << "Error while parsing command line options:" << endl
      << "The required option " << e.get_option_name() << " was not specified." << endl
      << default_error_msg << endl;
    return(-1);
  }

  if(config_path != "") {
    ifstream ifs(config_path.c_str());
    if(!ifs) {
      cout << "Error: can not open config file: " << config_path << endl
        << default_error_msg << endl;
      return(-1);
    } else {
      try{
        store(parse_config_file(ifs, config_file_options), vm);
      } catch(po::multiple_occurrences e) {
        cout << "Error while parsing config file:" << endl
          << "Option " << e.get_option_name() << " was specified multiple times." << endl
          << default_error_msg << endl;
        return(-1);
      } catch(po::unknown_option e) {
        cout << "Error while parsing config file:" << endl
          << "Option " << e.get_option_name() << " not known." << endl
          << default_error_msg << endl;
        return(-1);
      } catch(po::invalid_option_value e) {
        cout << "Error while parsing config file:" << endl
          << "The value specified for option " << e.get_option_name() << " has an invalid format." << endl
          << default_error_msg << endl;
        return(-1);
      }
      notify(vm);
    }
  }

  // generate an output path stem if the user did not specify one
  if(not vm.count("output")) {
    options.label = generate_random_label(options.exec_info.program_name, 0, options.verbosity);
    while(boost::filesystem::exists(options.label + ".hmm"))
      options.label = generate_random_label(options.exec_info.program_name, 5, options.verbosity);
    if(options.verbosity >= Verbosity::info)
      cout << "Using \"" << options.label << "\" as label to generate output file names." << endl;
  }

  // set the number of threads with OpenMP
  if(vm.count("threads"))
    omp_set_num_threads(options.n_threads);

  // print information about specified motifs, paths, and objectives
  if(options.verbosity >= Verbosity::debug) {
    cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
    cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
    cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
  }

  // check and harmonize specified motifs, paths, and objectives
  Specification::harmonize(options.motif_specifications, options.paths, options.objectives);

  // print information about specified motifs, paths, and objectives
  if(options.verbosity >= Verbosity::debug) {
    cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
    cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
    cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
  }

  // ensure that the user specified a positive number, if anything, for sampling min_size
  if(vm.count("smin")) {
    if(options.sampling.min_size < 0) {
      cout << "Please note that the minimal length for MCMC sampling must be non-negative." << endl
        << default_error_msg << endl;
      return(-1);
    }
  } else
    options.sampling.min_size = -1;

  // ensure that the user specified a positive number, if anything, for sampling max_size
  if(vm.count("smax")) {
    if(options.sampling.max_size < 0) {
      cout << "Please note that the maximal length for MCMC sampling must be non-negative." << endl
        << default_error_msg << endl;
      return(-1);
    } else if(options.sampling.max_size < options.sampling.min_size) {
      cout << "Please note that the maximal length for MCMC sampling must be larger than the minimal size." << endl
        << default_error_msg << endl;
      return(-1);
    }
  } else
    options.sampling.max_size = -1;

  // Initialize the plasma options
  fixup_seeding_options(options);

  // deprecation warning for option --hmmscore
  if(hmm_score_seeding) {
    options.seeding.keep_all = true;
    std::cout << "Warning: option --hmmscore is deprecated. It has been superseded by the --keepall option." << std::endl;
  }

  if(options.seeding.keep_all)
    options.model_choice = ModelChoice::HMMScore;
  else
    options.model_choice = ModelChoice::SeedScore;
  if(options.termination.past == 0) {
    std::cout << "Error: the value of --past must be a number greater than 0." << std::endl;
  }

  if(options.simultaneity == Training::Simultaneity::Simultaneous and options.wiggle != 0 and options.model_choice != ModelChoice::HMMScore) {
    std::cout << "Warning: option --wiggle deactivated due to simultaneous motif seeding based on seed scores." << std::endl;
    std::cout << "If you want to use wiggle variants please use --hmmscore." << std::endl;
    options.wiggle = 0;
  }

  if(options.termination.max_iter == 0 and options.sampling.do_sampling) {
    options.termination.max_iter = 1000;
    std::cout << "Note: did not specify the number of iterations to perform (--maxiter)." << std::endl
      << "We will now do " << options.termination.max_iter << " iterations." << std::endl;
  }

  bool any_named = false;
  for(auto &x: options.paths)
    if(not x.motifs.empty()) {
      any_named = true;
      break;
    }
  if(not any_named)
    for(auto &x: options.paths)
      for(auto &s: options.motif_specifications)
        x.motifs.insert(s.name);

  if(options.load_paths.empty() and options.motif_specifications.empty()) {
    std::cout << "Error: you must either specify at least one of:" << std::endl
      << "1. a path from which to load HMM parameter (--load)" << std::endl
      << "2. one or more motifs (--motif)" << std::endl
      << default_error_msg << endl;
    return(-1);
  }

  if(not options.long_names)
    if(options.seeding.keep_all or options.wiggle > 0) {
      std::cout << "Warning: you specified HMM score seed selection and / or wiggle variants, but did not specify --longnames." << std::endl
        << "Adding option --longnames." << std::endl;
      options.long_names = true;
    }

  if(options.line_search.eta <= options.line_search.mu) {
    std::cout << "Error: the Moré-Thuente η parameter must be larger than the µ parameter." << std::endl;
    exit(-1);
  }


  // initialize RNG
  if(options.verbosity >= Verbosity::info)
    std::cout << "Initializing random number generator with salt " << options.random_salt << "." << std::endl;
  srand(options.random_salt);

  // main routine
  perform_analysis(options);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  if(options.verbosity >= Verbosity::info) {
    std::cout << "CPU time used = " << cpu_time_used << std::endl;
    if(options.verbosity >= Verbosity::verbose) {
      tms tm;
      clock_t some_time = times(&tm);
      std::cout << "times() return " << ((double) some_time) << std::endl;
      std::cout << "utime = " << tm.tms_utime << std::endl;
      std::cout << "stime = " << tm.tms_stime << std::endl;
      std::cout << "cutime = " << tm.tms_cutime << std::endl;
      std::cout << "cstime = " << tm.tms_cstime << std::endl;
    }
  }

  return(EXIT_SUCCESS);
}

