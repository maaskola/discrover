/*
 * =====================================================================================
 *
 *       Filename:  cli.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include "options.hpp"
#include "../random_seed.hpp"
#include "../logo/cli.hpp"
#include "../executioninformation.hpp"

using namespace std;

string form_switch(const string &prefix, const string &s,
                   bool add_short = true) {
  if (prefix == "") {
    if (add_short)
      return s + "," + s[0];
    else
      return s;
  } else
    return prefix + s;
}

boost::program_options::options_description gen_plasma_options_description(
    Seeding::Options &options, size_t cols, const string &prefix,
    const string &name, bool include_all, bool allow_short) {
  namespace po = boost::program_options;
  po::options_description desc(name, cols);

  Seeding::Objective default_objective;
  default_objective.measure = Measures::Discrete::Measure::MutualInformation;

  if (include_all)
    desc.add_options()
      ("fasta,f", po::value(&options.paths)->required(),
       "FASTA file definition consisting of up to three colon-separated terms. "
       "May be given multiple times. "
       "Syntax as follows:\n"
       "\n"
       "[NAMES:[CONTRAST:]]PATH\n"
       "\n"
       "NAMES    \tcomma-separated list of names of motifs enriched in this file, or 'control'; directional objective functions (see description of score specifications) need to know which of the FASTA files are to treated as signal sequences; thus for these it is important to name the motifs and specify them in the FASTA file specification\n"
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
       // TODO note usage of the 'control' motif name
      )
      ("motif,m", po::value(&options.motif_specifications),
       "Motif definition consisting of up to two colon-separated terms. "
       "May be given multiple times. "
       "Syntax as follows:\n"
       "\n"
       "[NAME:]MOTIFSPEC\n"
       "\n"
       "NAME      \tunique name for the motif; may not be 'control'; directional objective functions (see description of score specifications) need to know which of the FASTA files are to treated as signal sequences; thus for these it is important to name the motifs and specify them in the FASTA file specification\n"
       "MOTIFSPEC \tmotif specification; see below\n"
       "\n"
       "NAME is optional. "
       "If NAME is not given, a name of the form 'motifN' is generated, where N are numbers increasing from 0.\n"
       "\n"
       "MOTIFSPEC can be given in multiple ways:\n"
       "1. \tUsing the IUPAC code for nucleic acids.\n"
       "2. \tA length specification. "
       "Motifs of the indicated lengths are sought. "
       "A length specification is a comma separated list of length ranges, "
       "where a length range is either a single length, or an expression of the form 'M-N' to indicate lengths M up to N. "
       "A length specification also allows to specify a multiplicity separated by an 'x'. "
       "Thus examples are '8' for just length 8, '5-7' for lengths 5, 6, and 7, '5-8x2' for two motifs of lengths 5 to 8, '5-8,10x3' for three motifs of lengths 5, 6, 7, 8, and 10."
       "\n"
       )
       ("score", po::value(&options.objectives)->default_value(Seeding::Objectives(1,default_objective), "mi"),
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
        "signal_freq  \tRelative occurrence frequency in the signal sequences\n"
        "control_freq \tRelative occurrence frequency in the control sequences\n"
        "mi           \tMutual information of condition and motif occurrence (MICO)\n"
        // "ri           \tRank mutual information\n"
        "mcc          \tMatthews correlation coefficient; directional measure (see note below)\n"
        "dfreq        \tDifference of frequency of sequences with motif occurrences, similar to the motif discovery methods DIPS and DECOD, see doi: 10.1093/bioinformatics/btl227 and 10.1093/bioinformatics/btr412; directional measure (see note below)\n"
        "gtest        \tG-test for association\n"
        "gtest_logp   \tlog-P value of G-test for association, corrected for multiple testing\n"
        // "gtest_logp_raw \tlog-P value of G-test for association, not corrected for multiple testing\n"
        "\n"
        "Note: contrasts for score specifications that use discriminative significance measures have to comprise at least two FASTA file specifications, otherwise these contrasts will be extended into binary contrasts by generating sequence shuffles to be used as control sequences.\n"
        "\n"
        "Note that the directional measures (mcc and dfreq) require you to specify in which of the samples the motif is expected. "
        "That is to say: you have to give a name to your motif in the motif specification, and add that name to the list of names given in a FASTA file specification."
        "\n"
        )
      ("revcomp,r", po::bool_switch(&options.revcomp), "Respect motif occurrences on the reverse complementary strand. Useful for DNA sequence motif analysis. Default is to consider only occurrence on the forward strand.")
      ("nseq", po::value(&options.n_seq)->default_value(0), "Use only the first N sequences of each file. Use 0 to indicate all sequences.")
      ;
  else {
    options.paths = {};
    options.motif_specifications = {};
    options.objectives = {1, {"mi"}};
    options.revcomp = false;
    options.n_seq = 0;
  }

  if (not include_all)
    desc.add_options()
      (form_switch(prefix, "seedscore", allow_short).c_str(), po::value(&options.objectives)->default_value(Seeding::Objectives(1,Seeding::Objective("mi")), "mi"), "Which objective function to evaluate. Available are 'signal_freq', 'control_freq', 'mi', 'mcc', 'dfreq', 'gtest', 'gtest_logp', 'gtest_logp_raw'.")
      ;

  desc.add_options()
    (form_switch(prefix, "algo", false).c_str(), po::value(&options.algorithm)->default_value(Seeding::Algorithm::Plasma, "plasma"), (string() + "Seeding algorithm. Available are 'plasma', 'mcmc', " + (DREME_FOUND ? "'dreme', " : "") + "and 'all'. Multiple algorithms can be used by separating them by comma.").c_str())
    (form_switch(prefix, "any", false).c_str(), po::bool_switch(&options.no_enrichment_filter), "Whether to allow motifs enriched in the opposite direction.")
    (form_switch(prefix, "filter", false).c_str(), po::value(&options.occurrence_filter)->default_value(Seeding::OccurrenceFilter::MaskOccurrences, "mask"), "How to filter motif occurrences upon identifying a motif. Available are 'remove' and 'mask'.")
    (form_switch(prefix, "cand", false).c_str(), po::value(&options.plasma.max_candidates)->default_value(100), "How many candidates to maintain.")
    (form_switch(prefix, "deg", allow_short).c_str(), po::value(&options.plasma.degeneracies), "Which degrees of degeneracy to consider. May be given multiple times. A sequence of length N has a maximal degeneracy of 3*N. Unlimited if unspecified.")
    (form_switch(prefix, "rdeg", false).c_str(), po::value(&options.plasma.rel_degeneracy)->default_value(1), "Limit relative degeneracy. 1 corresponds to full degeneracy, and 0 to no degeneracy. For a sequence of length N the degeneracy is maximally 3*N. Thus for a motif of length 8 a maximal relative degeneracy of 0.2 allows (rounded down) 4 degrees of degeneracy.")
    (form_switch(prefix, "generalize", allow_short).c_str(), po::bool_switch(&options.plasma.per_degeneracy), "Whether to report the best motifs at each level of degeneracy. Default is to report only the best motif across all levels of degeneracy. In addition, the best motifs of levels of degeneracy given by --deg will be reported.")
    (form_switch(prefix, "best", false).c_str(), po::bool_switch(&options.only_best), "Whether to report only the single best motif for each motif specification.  Default is to report for each motif specification the best result for each length.")
    (form_switch(prefix, "strict", false).c_str(), po::bool_switch(&options.strict), "Do not allow insignificant seeds.")
    (form_switch(prefix, "fix_mspace", false).c_str(), po::bool_switch(&options.fixed_motif_space_mode), "Deactivate dynamic motif space mode. Influences how the multiple-testing correction for the log-p value of the G-test is calculated.")
    (form_switch(prefix, "allowIUPAC", false).c_str(), po::bool_switch(&options.allow_iupac_wildcards), "Interpret IUPAC wildcard symbols in FASTA files. When this option is used e.g. S (strong) matches C and G, and so on. Importantly, N matches any character! Use non-IUPAC characters for positions where the sequence is unknown or masked, e.g. you could use '-' for this. By default, only A, C, G, and T characters (and their lower case variants) are encoded while all other characters are interpreted as masked.")
    ;
  if (include_all)
    desc.add_options()
      ("weight", po::bool_switch(&options.weighting), "When combining objective functions across multiple contrasts, combine values by weighting with the number of sequences per contrasts.")
      ("pcount", po::value(&options.pseudo_count)->default_value(1), "The number of pseudo counts to add to each cell of contingency tables.")
      ("word,w", po::bool_switch(&options.word_stats), "Perform nucleotide level statistics instead of on sequence level.")
      ("time", po::bool_switch(&options.measure_runtime), "Output information about how long certain parts take to execute.")
      ("print", po::bool_switch(&options.dump_viterbi), "Print out sequences annotated with motif occurrences.")
      ("bed", po::bool_switch(&options.dump_bed), "Generate a BED file with positions of motif occurrence.")
      ("threads", po::value(&options.n_threads), "Number of threads. If not given, as many are used as there are CPU cores on this machine.")
      ("output,o", po::value(&options.label),
       "Output file names are generated from this label. If not given, the output label will be 'plasma_XXX' where XXX is a string to make the label unique. The output files comprise:\n"
       "If --pdf or -png are used, sequence logos of the found motifs are generated with file names based on this output label.")
      ("salt", po::value(&options.mcmc.random_salt), "Seed for the pseudo random number generator (used e.g. for sequence shuffle generation and MCMC sampling). Set this to get reproducible results.")
      ;
  else {
    options.pseudo_count = 0;
    options.word_stats = false;
    options.measure_runtime = false;
    options.dump_viterbi = false;
    options.dump_bed = false;
    options.label = generate_random_label("plasma", 0, options.verbosity);
    options.mcmc.random_salt = generate_rng_seed();
  }

  if (include_all) {
    po::options_description mcmc_desc("MCMC optimization options", cols);
    string mcmc_prefix = "mcmc_";
    mcmc_desc.add_options()
      ("partemp", po::value(&options.mcmc.n_parallel)->default_value(6), "Parallel chains to run for parallel tempering.")
      ("maxiter", po::value(&options.mcmc.max_iter)->default_value(1000), "Maximal number of iterations to perform during MCMC seeding.")
      ("temp", po::value(&options.mcmc.temperature)->default_value(1e-3), "When performing MCMC sampling use this temperature. The temperatures of parallel chains is decreasing by factors of two.")
      ;
    desc.add(mcmc_desc);

#if CAIRO_FOUND
    po::options_description logo_options
        = gen_logo_options_description(options.logo, Logo::CLI::IUPAC, cols);
    desc.add(logo_options);
#endif
  } else {
    options.mcmc.temperature = 1e-3;
    options.mcmc.n_parallel = 6;
    options.mcmc.max_iter = 1000;
  }

  return desc;
}
