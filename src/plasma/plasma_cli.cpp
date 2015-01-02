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

#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include "options.hpp"
#include "../random_seed.hpp"
#include "../executioninformation.hpp"

using namespace std;

string form_switch(const string &prefix, const string& s, bool add_short=true) {
  if(prefix == "") {
    if(add_short)
      return(s + "," + s[0]);
    else
      return(s);
  } else
    return(prefix + s);
}

boost::program_options::options_description gen_plasma_options_description(Seeding::Options &options,
    const string &prefix,
    const string &name,
    size_t cols,
    bool include_all,
    bool allow_short) {
  namespace po = boost::program_options;
  po::options_description desc(name, cols);
  if(include_all)
    desc.add_options()
      ("fasta,f", po::value<vector<Specification::Set>>(&options.paths)->required(), "FASTA file(s) with nucleic acid sequences.")
      ("motif,m", po::value<vector<Specification::Motif>>(&options.motif_specifications)->required(), "Motif specification. "
       "May be given multiple times, and can be specified in multiple ways:\n"
       "1. \tUsing the IUPAC code for nucleic acids.\n"
       "2. \tA length specification. "
       "Motifs of the indicated lengths are sought. "
       "A length specification is a comma separated list of length ranges, "
       "where a length range is either a single length, or an expression of the form 'x-y' to indicate lengths x up to y. "
       "A length specification also allows to specify a multiplicity separated by an 'x'. "
       "Thus examples are '8' for just length 8, '5-7' for lengths 5, 6, and 7, '5-8x2' for two motifs of lengths 5 to 8, '5-8,10x3' for three motifs of lengths 5, 6, 7, 8, and 10.\n"
       "Regardless of the way the motif is specified, it may be given a name. The syntax is [NAME:]MOTIFSPEC.")
      ("score,s", po::value<Seeding::Objectives>(&options.objectives)->default_value(Seeding::Objectives(1,Seeding::Objective("mi")), "mi"), "Which objective function to evaluate. TODO: documentation needs updating to reflect more advanced options for this argument. Available are 'signal_freq', 'control_freq', 'mi', 'mcc', 'delta_freq', 'gtest', 'gtest_logp', 'gtest_logp_raw'.")
      ("revcomp,r", po::bool_switch(&options.revcomp), "Also consider the reverse complements of the sequences.")
      ("nseq", po::value<size_t>(&options.n_seq)->default_value(0), "Use only the first N sequences of each file. Use 0 to indicate all sequences.")
      ;
  else {
    options.paths = {};
    options.motif_specifications = {};
    options.objectives = {1, {"mi"}};
    options.revcomp = false;
    options.n_seq = 0;
  }

  if(not include_all)
    desc.add_options()
      (form_switch(prefix, "seedscore", allow_short).c_str(), po::value<Seeding::Objectives>(&options.objectives)->default_value(Seeding::Objectives(1,Seeding::Objective("mi")), "mi"), "Which objective function to evaluate. TODO: documentation needs updating to reflect more advanced options for this argument. Available are 'signal_freq', 'control_freq', 'mi', 'mcc', 'delta_freq', 'gtest', 'gtest_logp', 'gtest_logp_raw'.")
      ;

  desc.add_options()
    (form_switch(prefix, "algo", allow_short).c_str(), po::value<Seeding::Algorithm>(&options.algorithm)->default_value(Seeding::Algorithm::Plasma, "plasma"), "Which algorithm to use for seeding. Available are 'plasma', 'mcmc', and 'all'. Multiple algorithms can be used by separating them by comma.")
    (form_switch(prefix, "any", false).c_str(), po::bool_switch(&options.no_enrichment_filter), "Whether to allow motifs enriched in the opposite direction.")
    (form_switch(prefix, "filter", false).c_str(), po::value<Seeding::OccurrenceFilter>(&options.occurrence_filter)->default_value(Seeding::OccurrenceFilter::MaskOccurrences, "mask"), "How to filter motif occurrences upon identifying a motif. Available are 'remove' and 'mask'.")
    (form_switch(prefix, "cand", allow_short).c_str(), po::value<size_t>(&options.plasma.max_candidates)->default_value(100), "How many candidates to maintain.")
    (form_switch(prefix, "deg", allow_short).c_str(), po::value<vector<size_t>>(&options.plasma.degeneracies), "Which degrees of degeneracy to consider. May be given multiple times. A sequence of length N has a maximal degeneracy of 3*N. Unlimited if unspecified.")
    (form_switch(prefix, "rdeg", false).c_str(), po::value<double>(&options.plasma.rel_degeneracy)->default_value(1), "Limit relative degeneracy. 1 corresponds to full degeneracy, and 0 to no degeneracy. For a sequence of length N the degeneracy is maximally 3*N. Thus for a motif of length 8 a maximal relative degeneracy of 0.2 allows (rounded down) 4 degrees of degeneracy.")
    (form_switch(prefix, "generalize", allow_short).c_str(), po::bool_switch(&options.plasma.per_degeneracy), "Whether to report the best motifs at each level of degeneracy. Default is to report only the best motif across all levels of degeneracy. In addition, the best motifs of levels of degeneracy given by --deg will be reported.")
    (form_switch(prefix, "best", false).c_str(), po::bool_switch(&options.only_best), "Whether to report only the single best motif for each motif specification.  Default is to report for each motif specification the best result for each length.")
    (form_switch(prefix, "strict", false).c_str(), po::bool_switch(&options.strict), "Do not allow insignificant seeds.")
    (form_switch(prefix, "fix_mspace", false).c_str(), po::bool_switch(&options.fixed_motif_space_mode), "Deactivate dynamic motif space mode. Influences how the multiple-testing correction for the log-p value of the G-test is calculated.")
    (form_switch(prefix, "allowIUPAC", false).c_str(), po::bool_switch(&options.allow_iupac_wildcards), "Interpret IUPAC wildcard symbols in FASTA files. When this option is used e.g. S (strong) matches C and G, and so on. Importantly, N matches any character! Use non-IUPAC characters for positions where the sequence is unknown or masked, e.g. you could use '-' for this. By default, only A, C, G, and T characters (and their lower case variants) are encoded while all other characters are interpreted as masked.")
    ;
  if(include_all)
    desc.add_options()
      ("weight", po::bool_switch(&options.weighting), "When combining objective functions across multiple contrasts, combine values by weighting with the number of sequences per contrasts.")
      ("pcount,p", po::value<double>(&options.pseudo_count)->default_value(1), "The number of pseudo counts to add to each cell of contingency tables.")
      ("word,w", po::bool_switch(&options.word_stats), "Perform nucleotide level statistics instead of on sequence level.")
      ("time,t", po::bool_switch(&options.measure_runtime), "Report running times.")
      ("print", po::bool_switch(&options.dump_viterbi), "Print out sequences annotated with motif occurrences.")
#if CAIRO_FOUND
      ("pdf", po::bool_switch(&options.pdf_logo), "Generate PDF files with sequence logos of the found motifs.")
      ("png", po::bool_switch(&options.png_logo), "Generate PNG files with sequence logos of the found motifs.")
#endif
      ("threads,T", po::value<size_t>(&options.n_threads)->default_value(omp_get_num_procs()), "The number of threads to use. If this in not specified, the value of the environment variable OMP_NUM_THREADS is used if that is defined, otherwise it will use as many as there are CPU cores on this machine.")
      ("output,o", po::value<string>(&options.label), "Output file names are generated from this label. If this option is not specified the output label will be 'plasma_XXX' where XXX is a string to make the label unique.")
      ("salt", po::value<unsigned int>(&options.mcmc.random_salt), "Seed for the pseudo random number generator (used e.g. for sequence shuffle generation and MCMC sampling). Set this to get reproducible results.")
      ;
  else {
    options.pseudo_count = 0;
    options.word_stats = false;
    options.measure_runtime = false;
    options.dump_viterbi = false;
    options.label = generate_random_label("plasma", 0, options.verbosity);
    options.mcmc.random_salt = generate_rng_seed();
#if CAIRO_FOUND
    options.pdf_logo = false;
    options.png_logo = false;
#endif
  }

  if(include_all) {
    po::options_description mcmc_desc("MCMC seeding options", cols);
    string mcmc_prefix = "mcmc_";
    mcmc_desc.add_options()
      ("temp", po::value<double>(&options.mcmc.temperature)->default_value(1e-3), "When performing Gibbs sampling use this temperature. The temperatures of parallel chains is decreasing by factors of two.")
      ("maxiter", po::value<size_t>(&options.mcmc.max_iter)->default_value(1000), "Maximal number of iterations to perform during MCMC seeding.")
      ("partemp", po::value<size_t>(&options.mcmc.n_parallel)->default_value(6), "Parallel chains to run for parallel tempering.")
      ;
    desc.add(mcmc_desc);
  } else {
    options.mcmc.temperature = 1e-3;
    options.mcmc.n_parallel = 6;
    options.mcmc.max_iter = 1000;
  }


  return(desc);
}

