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

#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include "options.hpp"

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

boost::program_options::options_description gen_iupac_options_description(Seeding::options_t &options,
    const string &prefix,
    const string &name,
    size_t cols,
    bool include_all,
    bool allow_short) {
  namespace po = boost::program_options;
  po::options_description desc(name, cols);
  if(include_all)
    desc.add_options()
      ("fasta,f", po::value<vector<Specification::DataSet>>(&options.paths)->required(), "FASTA file(s) with nucleic acid sequences.")
      ("motif,m", po::value<vector<Specification::Motif>>(&options.motif_specifications)->required(), "Motif specification. May be given multiple times, and can be specified in multiple ways:\n"
       "1. \tUsing the IUPAC code for nucleic acids.\n"
       "2. \tA length specification. Motifs of the indicated lengths are sought. A length specification is a comma separated list of length ranges, where a length range is either a single length, or an expression of the form 'x-y' to indicate lengths x up to y.\n"
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
    (form_switch(prefix, "nmotif", allow_short).c_str(), po::value<size_t>(&options.n_motifs)->default_value(1), "How many motifs to determine.")
    (form_switch(prefix, "any", allow_short).c_str(), po::bool_switch(&options.no_enrichment_filter), "Whether to allow motifs enriched in the opposite direction.")
    (form_switch(prefix, "filter", false).c_str(), po::value<Seeding::OccurrenceFilter>(&options.occurrence_filter)->default_value(Seeding::OccurrenceFilter::mask_occurrence, "mask"), "How to filter motif occurrences upon identifying a motif. Available are 'remove' and 'mask'.")
    (form_switch(prefix, "cand", allow_short).c_str(), po::value<size_t>(&options.max_candidates)->default_value(100), "How many candidates to maintain.")
    (form_switch(prefix, "deg", allow_short).c_str(), po::value<vector<size_t>>(&options.degeneracies), "Which degrees of degeneracy to consider. May be given multiple times. A sequence of length N has a maximal degeneracy of 3*N. Unlimited if unspecified.")
    (form_switch(prefix, "rdeg", false).c_str(), po::value<double>(&options.rel_degeneracy)->default_value(1), "Limit relative degeneracy. 1 corresponds to full degeneracy, and 0 to no degeneracy. For a sequence of length N the degeneracy is maximally 3*N. Thus for a motif of length 8 a maximal relative degeneracy of 0.2 allows (rounded down) 4 degrees of degeneracy.")
    (form_switch(prefix, "generalize", allow_short).c_str(), po::bool_switch(&options.per_degeneracy), "Whether to report the best motifs at each level of degeneracy. Default is to report only the best motif across all levels of degeneracy. In addition, the best motifs of levels of degeneracy given by --deg will be reported.")
    (form_switch(prefix, "keepall", false).c_str(), po::bool_switch(&options.keep_all), "Whether to report for each motif specification the best result for each length. Default is to report only the single best motif for each motif specification.")
    (form_switch(prefix, "strict", false).c_str(), po::bool_switch(&options.strict), "Do not allow insignificant seeds.")
    ;
  if(include_all)
    desc.add_options()
      ("pcount,p", po::value<double>(&options.pseudo_count)->default_value(1), "The number of pseudo counts to add to each cell of contingency tables.")
      ("word,w", po::bool_switch(&options.word_stats), "Perform nucleotide level statistics instead of on sequence level.")
      ("time,t", po::bool_switch(&options.measure_runtime), "Report running times.")
      ("print", po::bool_switch(&options.dump_viterbi), "Print out sequences annotated with motif occurrences.")
      ("threads,T", po::value<size_t>(&options.n_threads)->default_value(omp_get_num_procs()), "The number of threads to use. If this in not specified, the value of the environment variable OMP_NUM_THREADS is used if that is defined, otherwise it will use as many as there are CPU cores on this machine.")
      ;
  else {
    options.pseudo_count = 0;
    options.word_stats = false;
    options.measure_runtime = false;
    options.dump_viterbi = false;
  }
  return(desc);
}

