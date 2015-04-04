/* =====================================================================================
 * Copyright (c) 2015, Jonas Maaskola
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
 *       Filename:  cli.cpp
 *
 *    Description:  Command line options generation for Discrover
 *
 *        Created:  Fri Jan 30 20:40:22 2015 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "../plasma/cli.hpp"
#include "../logo/cli.hpp"
#include "cli.hpp"

using namespace std;

namespace po = boost::program_options;

void gen_discrover_cli(size_t cols, string &config_path, Options::HMM &options,
                       po::options_description &common_options,
                       po::options_description &visible_options,
                       po::options_description &hidden_options,
                       po::options_description &cmdline_options,
                       po::options_description &config_file_options) {
  // Declare the supported options.
  po::options_description generic_options("Generic options", cols);
  po::options_description basic_options("Basic options, required", cols);
  po::options_description basic_options_optional("Basic options", cols);
  po::options_description init_options("Initialization options", cols);
  po::options_description eval_options("Evaluation options", cols);
  po::options_description advanced_options("Advanced options", cols);
  po::options_description conjugate_options("Conjugate gradient calculation options", cols);
  po::options_description multi_motif_options("Multiple motif mode options", cols);
  po::options_description mmie_options("MMIE options", cols);
  po::options_description linesearching_options("Line searching options", cols);
  po::options_description sampling_options("MCMC optimization options", cols);
  po::options_description termination_options("Termination options", cols);

  po::options_description seeding_options = gen_plasma_options_description(
      options.seeding, cols, "",
      "Seeding options for IUPAC regular expression finding", false,
      false);

#if CAIRO_FOUND
  po::options_description logo_options
      = gen_logo_options_description(options.logo, Logo::CLI::HMM, cols);
#endif

  generic_options.add_options()
    ("config", po::value(&config_path), "Read options from a configuration file. ")
    ("help,h", "Produce help message. Combine with -v or -V for additional commands.")
    ("version", "Print out the version. Also show git SHA1 with -v.")
    ("verbose,v", "Be verbose about the progress.")
    ("noisy,V", "Be very verbose about the progress.")
    ;

  basic_options.add_options()
    ("fasta,f", po::value(&options.paths)->required(),
     "FASTA file definition consisting of up to three colon-separated terms. "
     "May be given multiple times. "
     "Syntax as follows:\n"
     "\n"
     "[NAMES:[CONTRAST:]]PATH\n"
     "\n"
     "NAMES    \tcomma-separated list of names of motifs enriched in this file, or 'control'\n"
     "CONTRAST \tname of a contrast this file belongs to\n"
     "PATH     \tpath of FASTA file with sequences\n"
     "\n"
     "NAMES and CONTRAST are optional. "
     "All FASTA files for which no CONTRAST is given are associated with a common, automatically named contrast.\n"
     "\n"
     "FASTA files can be grouped into contrasts on which the --score definitions (see below) can act.\n"
     "\n"
     "If PATH contains at least one colon, please prepend colons to disambiguate.\n"
     "\n"
     "Note: usage of -f / --fasta is optional; all free arguments are taken to be paths of FASTA files."
     "\n"
     "FASTA files tagged 'control' are treated slightly differently: "
     "they are excluded from the set of sequences on which parameter re-estimation is performed (motif occurrence priors and background emissions will thus not be learned from these).\n"
     "Note: shuffled sequences are implicitly tagged 'control'.\n"
     "\n"
     // TODO note usage of the 'control' motif name
     )
    ("motif,m", po::value(&options.motif_specifications),
     "Motif definition consisting of up to three colon-separated terms. "
     "May be given multiple times. "
     "Syntax as follows:\n"
     "\n"
     "[NAME:[INSERT:]]MOTIFSPEC\n"
     "\n"
     "NAME      \tunique name for the motif; may not be 'control'\n"
     "INSERT    \tpositions after which insertions are allowed\n"
     "MOTIFSPEC \tmotif specification; see below\n"
     "\n"
     "NAME and INSERT are optional. "
     "If NAME is not given, a name of the form 'motifN' is generated, where N are numbers increasing from 0.\n"
     "\n"
     "The INSERT string is a comma-separated list of positions after which to add an insertion state. "
     "1-indexed; 1 ≤ position < motif length must hold for each position.\n"
     "\n"
     "MOTIFSPEC can be given in multiple ways:\n"
     "1. \tUsing the IUPAC code for nucleic acids.\n"
     "2. \tA length specification. "
     "Motifs of the indicated lengths are sought. "
     "A length specification is a comma separated list of length ranges, "
     "where a length range is either a single length, or an expression of the form 'M-N' to indicate lengths M up to N. "
     "A length specification also allows to specify a multiplicity separated by an 'x'. "
     "Thus examples are '8' for just length 8, '5-7' for lengths 5, 6, and 7, '5-8x2' for two motifs of lengths 5 to 8, '5-8,10x3' for three motifs of lengths 5, 6, 7, 8, and 10.\n"
     "3. \tBy specifying the path to a file with a position weight matrix (PWM). "
     "If the path contains at least one colon, please prepend colons to disambiguate."
     "\n"
     )
   ;

  Training::Objective default_objective;
  default_objective.measure = Measure::MutualInformation;

  basic_options_optional.add_options()
    ("score", po::value(&options.objectives)->default_value(Training::Objectives(1,default_objective), "mi"),
     "Score definition consisting of up to three colon-separated terms. "
     "May be given multiple times. "
     "Syntax as follows:\n"
     "\n"
     "[MOTIF:[CONTRASTS:]]MEASURE\n"
     "\n"
     "MOTIF     \tname of the motif to optimize\n"
     "CONTRASTS \tcontrast specification\n"
     "MEASURE   \tsignificance measure to use\n"
     "See below for details.\n"
     "\n"
     "MOTIF and CONTRASTS are optional.\n"
     "\n"
     "If no motif name is given, all motifs (that are not named in any other score specification) are optimized with this score specification. "
     "It is an error when more than one score specification is given for any motif.\n"
     "\n"
     "The contrast specification is a set of contrast names prefixed by '+' or '-' to indicate whether the score on the contrast is to count positively or negatively. "
     "The first contrast name, if lacking a '+' or '-', is taken to count positively. "
     "If no contrast specification is given, the score definition uses all contrasts, counting each one positively.\n"
     "\n"
     "The following significance measures are available:\n"
     "none    \tDo not perform learning\n"
     "bw      \tLikelihood using Baum-Welch\n"
     "viterbi \tViterbi learning\n"
     "mi      \tMutual information of condition and motif occurrence (MICO)\n"
     "ri      \tRank mutual information\n"
     "mmie    \tMaximum mutual information estimation (MMIE), a.k.a. posterior classification probability\n"
     "mcc     \tMatthews correlation coefficient\n"
     "dlogl   \tLog-likelihood difference, like in the motif discovery method DME, see doi: 10.1073/pnas.0406123102\n"
     "dfreq   \tDifference of frequency of sequences with motif occurrences, similar to the motif discovery methods DIPS and DECOD, see doi: 10.1093/bioinformatics/btl227 and 10.1093/bioinformatics/btr412\n"
     "\n"
     "Note: contrasts for score specifications that use discriminative significance measures have to comprise at least two FASTA file specifications, otherwise these contrasts will be extended into binary contrasts by generating sequence shuffles to be used as control sequences."
     "\n"
     )
    ("revcomp,r", po::bool_switch(&options.revcomp), "Respect motif occurrences on the reverse complementary strand. Useful for DNA sequence motif analysis. Default is to consider only occurrence on the forward strand.")
    ;

  advanced_options.add_options()
    ("output,o", po::value(&options.label), (
      "Output file names are generated from this label. If not given, the output label will be '" + options.exec_info.program_name + "_XXX' where XXX is a string to make the label unique. The output files comprise:\n"
      ".hmm     \tParameter of the trained HMM. May be loaded later with --learn.\n"
      ".summary \tSummary information with number of occurrences of the motifs in each sample, and various generative and discriminative statistics.\n"
      ".viterbi \tFASTA sequences annotated with the Viterbi path, and sequence level statistics.\n"
      ".bed     \tBED file of coordinates of matches to the motifs in all sequences.\n"
      ".table   \tCoordinates and sequences of matches to the motifs in all sequences (extends the .bed output file).\n"
      ".pdf     \tSequence logos of the found motifs in PDF format.\n"
      ".png     \tSequence logos of the found motifs in PNG format.\n"
      "Note that, depending on the argument of --compress, the .viterbi, .bed, and .table files may be compressed, and require decompression for inspection.\n"
     ).c_str())
    ("threads", po::value(&options.n_threads), "Number of threads. If not given, as many are used as there are CPU cores on this machine.")
    ("time", po::bool_switch(&options.timing_information), "Output information about how long certain parts take to execute.")
    ("cv", po::value(&options.cross_validation_iterations)->default_value(0), "Number of cross validation iterations to do.")
    ("cv_freq", po::value(&options.cross_validation_freq)->default_value(0.9, "0.9"), "Fraction of data samples for training in cross validation.")
    ("nseq", po::value(&options.n_seq)->default_value(0), "Use only the first N sequences of each file. Use 0 to indicate all sequences.")
    ("iter", po::value(&options.termination.max_iter)->default_value(1000), "Maximal number of iterations to perform in training. A value of 0 means no limit, and that the training is only terminated by the tolerance.")
    ("salt", po::value(&options.random_salt), "Seed for the pseudo random number generator (used e.g. for sequence shuffle generation and MCMC sampling). Set this to get reproducible results.")
    ("weight", po::bool_switch(&options.weighting), "When combining objective functions across multiple contrasts, combine values by weighting with the number of sequences per contrasts.")
    ;

  conjugate_options.add_options()
    ("cg_mode", po::value(&options.conjugate.mode)->default_value(Options::Conjugate::Mode::None, "none"),
     "Conjugate gradient calculation method\n"
     "none \tNo conjugate gradient (steepest ascent)\n"
     "fr   \tFletcher-Reeves\n"
     "pr   \tPolak-Ribière\n"
     "hs   \tHestenes-Stiefel\n"
     "dy   \tDai-Yuan")
    ("cg_iter", po::value(&options.conjugate.restart_iteration)->default_value(0), "Number of iterations after which to reset conjugate gradient. Use 0 to never reset.")
    ("cg_thresh", po::value(&options.conjugate.restart_threshold)->default_value(0), "Threshold for gradient orthogonality (between 0 and 1) below which the conjugate gradient will be reset. Use 0 to never reset.")
    ;

  init_options.add_options()
    ("load,l", po::value(&options.load_paths), "Load HMM parameters from a .hmm file produced by an earlier run. Can be specified multiple times; then the first parameter file will be loaded, and motifs of the following parameter files are added.")
    ("selftrans", po::bool_switch(&options.self_transition), "Add self-transition edges to insert positions.")
    ("alpha", po::value(&options.alpha)->default_value(0.03, "0.03"), "Probability of alternative nucleotides. The nucleotides not included in the IUPAC character will have this probability.")
    ("lambda", po::value(&options.lambda)->default_value(1), "Initial value for prior with which a motif is expected.")
    ("wiggle", po::value(&options.wiggle)->default_value(0), "For automatically determined seeds, consider variants shifted up- and downstream by up to the specified number of positions.")
    ("extend", po::value(&options.extend)->default_value(0), "Extend seeds by this many Ns up- and downstream before HMM training.")
    ("padl", po::value(&options.left_padding)->default_value(0), "Add this many Ns upstream (to the left) of the seed.")
    ("padr", po::value(&options.right_padding)->default_value(0), "Add this many Ns downstream (to the right) of the seed.")
    ;

  eval_options.add_options()
    ("posterior", po::bool_switch(&options.evaluate.print_posterior), "During evaluation also print out the motif posterior probability.")
    ("condmotif", po::bool_switch(&options.evaluate.conditional_motif_probability), "During evaluation compute for every position the conditional motif likelihood considering only the motif emissions.")
    ("nosummary", po::bool_switch(&options.evaluate.skip_summary), "Do not print summary information.")
    ("noviterbi", po::bool_switch(&options.evaluate.skip_viterbi_path), "Do not print the Viterbi path.")
    ("nobed", po::bool_switch(&options.evaluate.skip_bed), "Do not generate BED file with positions of motif occurrences.")
    ("notable", po::bool_switch(&options.evaluate.skip_occurrence_table), "Do not print the occurrence table.")
    ("ric", po::bool_switch(&options.evaluate.perform_ric), "Perform a rank information coefficient analysis.")
    ;

  multi_motif_options.add_options()
    ("multiple", po::bool_switch(&options.multi_motif.accept_multiple), "Accept multiple motifs as long as the score increases. This can only be used with the objective function MICO.")
    ("relearn", po::value(&options.multi_motif.relearning)->default_value(Options::MultiMotif::Relearning::Full, "full"), "When accepting multiple motifs, whether and how to re-learn the model after a new motif is added. Choices: 'none', 'reest', 'full'.")
    ("resratio", po::value(&options.multi_motif.residual_ratio)->default_value(5.0), "Cutoff to discard new motifs in multi motif mode. The cutoff is applied on the ratio of conditional mutual information of the new motif and the conditions given the previous motifs. Must be non-negative. High values discard more motifs, and lead to less redundant motifs.")
    ;

  mmie_options.add_options()
    ("classp", po::value(&options.class_prior)->default_value(0.5),
     "Initial class prior.")
    ("motifp1", po::value(&options.conditional_motif_prior1)->default_value(0.6, "0.6"),
     "Initial conditional motif prior for the signal class.")
    ("motifp2", po::value(&options.conditional_motif_prior2)->default_value(0.03, "0.03"),
     "Initial conditional motif prior for the control class.")
    ("noclassp", po::bool_switch(&options.dont_learn_class_prior),
     "Don't learn the class prior.")
    ("nomotifp", po::bool_switch(&options.dont_learn_conditional_motif_prior),
     "Don't learn the conditional motif prior.")
    ;


  linesearching_options.add_options()
    ("LSmu", po::value(&options.line_search.mu)->default_value(0.1, "0.1"), "The parameter µ for the Moré-Thuente line search algorithm.")
    ("LSeta", po::value(&options.line_search.eta)->default_value(0.5, "0.5"), "The parameter η for the Moré-Thuente line search algorithm.")
    ("LSdelta", po::value(&options.line_search.delta)->default_value(0.66, "0.66"), "The parameter delta for the Moré-Thuente line search algorithm.")
    ("LSnum", po::value(&options.line_search.max_steps)->default_value(10, "10"), "How many gradient and function evaluation to perform maximally per line search.")
    ;

  termination_options.add_options()
    ("gamma", po::value(&options.termination.gamma_tolerance)->default_value(1e-4, "1e-4"), "Tolerance for the reestimation-type learning methods. Training stops when the L1 norm of the parameter change between iterations is less than this value.")
    ("delta", po::value(&options.termination.delta_tolerance)->default_value(1e-4, "1e-4"), "Relative score difference criterion tolerance for training algorithm termination: stops iterations when (f - f') / f < delta, where f' is the objective value of the past iteration, and f is the objective value of the current iteration.")
    ("epsilon", po::value(&options.termination.epsilon_tolerance)->default_value(0), "Gradient norm criterion tolerance for training algorithm termination: stops when ||g|| < epsilon * max(1, ||g||), where ||.|| denotes the Euclidean (L2) norm.")
    ("past", po::value(&options.termination.past)->default_value(1), "Distance for delta-based convergence test. This parameter determines the distance, in iterations, to compute the rate of decrease of the objective function.")
    ;

  sampling_options.add_options()
    ("sampling", po::bool_switch(&options.sampling.do_sampling), "Perform Monte-Carlo Markov chain (MCMC) sampling for parameter inference instead of re-estimation or gradient learning.")
    ("partemp", po::value(&options.sampling.n_parallel)->default_value(6), "Number of chains in parallel tempering.")
    ("temp", po::value(&options.sampling.temperature)->default_value(1e-3), "When performing MCMC sampling use this temperature. The temperatures of parallel chains is decreasing by factors of two.")
    ("smin", po::value(&options.sampling.min_size), "Minimal motif length for MCMC sampling. When unspecified defaults to initial motif length.")
    ("smax", po::value(&options.sampling.max_size), "Maximal motif length for MCMC sampling. When unspecified defaults to initial motif length.")
    ("nindel", po::value(&options.sampling.n_indels)->default_value(5), "Maximal number of positions that may be added or removed at a time. Adding and removing of nucleotides happens in front of and at the ends of the motif.")
    ("nshift", po::value(&options.sampling.n_shift)->default_value(5), "Maximal number of positions that the motif may be shifted by.")
    ;

  hidden_options.add_options()
    ("nosave", po::bool_switch(&options.dont_save_shuffle_sequences), "Do not save generated shuffle sequences.")
    ("bglearn", po::value(&options.bg_learning)->default_value(Training::Method::Reestimation, "em"), "How to learn the background. Available are 'fixed', 'em', 'gradient', where the 'em' uses re-estimation to maximize the likelihood contribution of the background parameters, while 'gradient' uses the discriminative objective function.")
    ("pscnt", po::value(&options.contingency_pseudo_count)->default_value(1.0, "1"), "The pseudo count to be added to the contingency tables in the discriminative algorithms.")
    ("pscntE", po::value(&options.emission_pseudo_count)->default_value(1.0, "1"), "The pseudo count to be added to the expected emission probabilities before normalization in the Baum-Welch algorithm.")
    ("pscntT", po::value(&options.transition_pseudo_count)->default_value(0.0, "0"), "The pseudo count to be added to the expected transition probabilities before normalization in the Baum-Welch algorithm.")
    ("compress", po::value(&options.output_compression)->default_value(Options::Compression::gzip, "gz"), "Compression method for larger output files. Available are: 'none', 'gz' or 'gzip', 'bz2' or 'bzip2'.") // TODO make the code conditional on the presence of zlib
    ("miseeding", po::bool_switch(&options.use_mi_to_seed), "Disregard automatic seeding choice and use MICO for seeding.")
    ("absthresh", po::bool_switch(&options.termination.absolute_improvement), "Whether improvement should be gauged by absolute value. Default is relative to the current score.")
    ("intermediate", po::bool_switch(&options.store_intermediate), "Write out intermediate parameters during training.")
    ("limitlogp", po::bool_switch(&options.limit_logp), "Do not report corrected log-P values greater 0 but report 0 in this case.")
    ;

  advanced_options
    .add(conjugate_options)
    .add(multi_motif_options)
    .add(mmie_options)
    .add(sampling_options);

  hidden_options
    .add(linesearching_options)
    .add(eval_options)
    .add(termination_options);

  po::options_description simple_options;
  simple_options
    .add(basic_options)
    .add(basic_options_optional);

  common_options
    .add(advanced_options)
    .add(seeding_options)
    .add(init_options);
#if CAIRO_FOUND
  common_options
    .add(logo_options);
#endif

  visible_options
    .add(generic_options)
    .add(simple_options);

  cmdline_options
    .add(generic_options)
    .add(simple_options)
    .add(common_options)
    .add(hidden_options);

  config_file_options
    .add(simple_options)
    .add(common_options)
    .add(hidden_options);
}
