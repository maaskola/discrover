/*
 * =====================================================================================
 *
 *       Filename:  options.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include <boost/algorithm/string.hpp>
#include "options.hpp"
#include "../random_seed.hpp"

using namespace std;

namespace Seeding {

Options::Options()
    : paths(),
      motif_specifications({}),
      objectives(),
      algorithm(Algorithm::Plasma),
      plasma(Plasma()),
      n_threads(1),
      revcomp(false),
      strict(false),
      pseudo_count(1),
      weighting(false),
      n_seq(0),
      word_stats(false),
      measure_runtime(false),
      occurrence_filter(OccurrenceFilter::RemoveSequences),
      only_best(false),
      verbosity(Verbosity::info),
      dump_viterbi(false),
      no_enrichment_filter(false),
      fixed_motif_space_mode(false),
      allow_iupac_wildcards(false),
      label(""){};

Options::Plasma::Plasma()
    : max_candidates(100),
      degeneracies(),
      rel_degeneracy(1),
      per_degeneracy(false){};

Options::MCMC::MCMC()
    : max_iter(1000),
      temperature(1e-3),
      n_parallel(6),
      random_salt(generate_rng_seed()){};

istream &operator>>(istream &in, OccurrenceFilter &filter) {
  string token;
  in >> token;
  if (token == "remove")
    filter = OccurrenceFilter::RemoveSequences;
  else if (token == "mask")
    filter = OccurrenceFilter::MaskOccurrences;
  else
    throw Exception::InvalidOccurrenceFilter(token);
  return in;
}

ostream &operator<<(ostream &os, const OccurrenceFilter &filter) {
  switch (filter) {
    case OccurrenceFilter::RemoveSequences:
      os << "remove";
      break;
    case OccurrenceFilter::MaskOccurrences:
      os << "mask";
      break;
  }
  return os;
}

Algorithm parse_algorithm(const string &token_) {
  string token(token_);
  boost::algorithm::to_lower(token);
  if (token == "plasma")
    return Algorithm::Plasma;
  else if (token == "dreme")
    return Algorithm::ExternalDREME;
  else if (token == "mcmc")
    return Algorithm::MCMC;
  else if (token == "all")
    return Algorithm::Plasma | Algorithm::ExternalDREME | Algorithm::MCMC;
  else
    throw Exception::InvalidAlgorithm(token);
}

istream &operator>>(istream &in, Algorithm &algorithm) {
  string algorithms;
  in >> algorithms;
  bool first = true;
  size_t pos;
  do {
    pos = algorithms.find(",");
    Algorithm algo = parse_algorithm(algorithms.substr(0, pos));
    if (first) {
      first = false;
      algorithm = algo;
    } else
      algorithm = algorithm | algo;
    algorithms = algorithms.substr(pos + 1);
  } while (pos != string::npos);
  Algorithm algo = parse_algorithm(algorithms);
  if (first)
    algorithm = algo;
  else
    algorithm = algorithm | algo;
  return in;
}
ostream &operator<<(ostream &os, const Algorithm &algorithm) {
  bool first = true;
  if ((algorithm & Algorithm::Plasma) == Algorithm::Plasma) {
    os << "plasma";
    first = false;
  }
  if ((algorithm & Algorithm::ExternalDREME) == Algorithm::ExternalDREME) {
    os << (first ? "" : ",") << "dreme";
    first = false;
  }
  if ((algorithm & Algorithm::MCMC) == Algorithm::MCMC)
    os << (first ? "" : ",") << "mcmc";
  return os;
}

ostream &operator<<(ostream &out, const Objective &objective) {
  out << Specification::to_string(objective);
  return out;
}

Objective objective_for_motif(const Objectives &objectives,
                              const Specification::Motif &motif) {
  for (auto objective : objectives)
    if (objective.motif_name == motif.name)
      return objective;
  throw Exception::NoMatchingObjectiveFound(motif.name);
}

namespace Exception {
InvalidOccurrenceFilter::InvalidOccurrenceFilter(const string &token_)
    : exception(), token(token_) {};
const char *InvalidOccurrenceFilter::what() const noexcept {
  string msg = "Error: invalid occurrence filter type '" + token + "'.";
  return msg.c_str();
}
InvalidAlgorithm::InvalidAlgorithm(const string &token_)
    : exception(), token(token_) {};
const char *InvalidAlgorithm::what() const noexcept {
  string msg = "Error: invalid seeding algorithm '" + token + "'.\n"
     + "Please use one of 'plasma', 'dreme', 'mcmc', or 'all'.\n"
     + "It is also possible to use multiple algorithms by separating them by comma.";
  return msg.c_str();
}
NoMatchingObjectiveFound::NoMatchingObjectiveFound(const string &motif_)
    : exception(), motif(motif_) {};
const char *NoMatchingObjectiveFound::what() const noexcept {
  string msg = "Error: no objective found for motif '" + motif + "'.";
  return msg.c_str();
}
}
}
