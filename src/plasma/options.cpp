/*
 * =====================================================================================
 *
 *       Filename:  options.cpp
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

#include <iostream>
#include <boost/algorithm/string.hpp>
#include "options.hpp"

using namespace std;


namespace Plasma {

  options_t::options_t() :
    paths(),
    motif_specifications({}),
    objectives(),
    algorithm(Algorithm::Plasma),
    n_threads(1),
    revcomp(false),
    strict(false),
    pseudo_count(1),
    n_seq(0),
    word_stats(false),
    measure_runtime(false),
    n_motifs(1),
    occurrence_filter(OccurrenceFilter::remove_seq),
    max_candidates(100),
    degeneracies(),
    rel_degeneracy(1),
    per_degeneracy(false),
    keep_all(false),
    verbosity(Verbosity::info),
    dump_viterbi(false),
    no_enrichment_filter(false)
    { };

  istream &operator>>(istream &in, OccurrenceFilter &filter) {
    string token;
    in >> token;
    if(token == "remove")
      filter = OccurrenceFilter::remove_seq;
    else if(token == "mask")
      filter = OccurrenceFilter::mask_occurrence;
    else {
      cout << "Couldn't parse occurrence filter type '" << token << "'." << endl;
      exit(-1);
    }
    return(in);
  }

  ostream &operator<<(ostream &os, const OccurrenceFilter &filter) {
    switch(filter) {
      case OccurrenceFilter::remove_seq:
        os << "remove";
        break;
      case OccurrenceFilter::mask_occurrence:
        os << "mask";
        break;
    }
    return(os);
  }

  Algorithm parse_algorithm(const string &token_) {
    string token(token_);
    boost::algorithm::to_lower(token);
    if(token == "plasma")
      return(Algorithm::Plasma);
    else if(token == "fire")
      return(Algorithm::FIRE);
    else if(token == "all")
      return(Algorithm::Plasma | Algorithm::FIRE);
    else {
      cout << "Seeding algorithm '" << token_ << "' unknown." << endl
        << "Please use one of 'plasma', 'fire', or 'all'." << endl
        << "It is also possible to use multiple algorithms by separating them by comma." << endl;
      exit(-1);
    }
  }

  istream &operator>>(istream &in, Algorithm &algorithm) {
    string algorithms;
    in >> algorithms;
    bool first = true;
    size_t pos;
    do {
      pos = algorithms.find(",");
      Algorithm algo = parse_algorithm(algorithms.substr(0, pos));
      if(first) {
        first = false;
        algorithm = algo;
      } else
        algorithm = algorithm | algo;
      algorithms = algorithms.substr(pos+1);
    } while(pos != string::npos);
    Algorithm algo = parse_algorithm(algorithms);
    if(first)
      algorithm = algo;
    else
      algorithm = algorithm | algo;
    cout << "Parsed algorithm: " << algorithm << endl;
    return(in);
  }
  ostream &operator<<(ostream &os, const Algorithm &algorithm) {
    bool first = true;
    if((algorithm & Algorithm::Plasma) == Algorithm::Plasma) {
      os << "plasma";
      first = false;
    }
    if((algorithm & Algorithm::FIRE) == Algorithm::FIRE)
      os << (first ? "" : ",") << "fire";
    return(os);
  }

  ostream &operator<<(ostream &out, const Objectives &objectives) {
    bool first = true;
    for(auto &objective: objectives) {
      if(first)
        first = false;
      else
        out << " ";
      out << Specification::to_string(objective);
    }
    return(out);
  }
}

