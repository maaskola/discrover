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
#include "options.hpp"

using namespace std;


namespace Seeding {

  options_t::options_t() :
    paths(),
    motif_specifications({}),
    objectives(),
    n_threads(1),
    revcomp(false),
    strict(false),
    pseudo_count(1),
    n_seq(0),
    word_stats(false),
    measure_runtime(false),
    n_motifs(1),
    occurrence_filter(OccurrenceFilter::RemoveSequences),
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
      filter = OccurrenceFilter::RemoveSequences;
    else if(token == "mask")
      filter = OccurrenceFilter::MaskOccurrences;
    else {
      cout << "Couldn't parse occurrence filter type '" << token << "'." << endl;
      exit(-1);
    }
    return(in);
  }

  ostream &operator<<(ostream &os, const OccurrenceFilter &filter) {
    switch(filter) {
      case OccurrenceFilter::RemoveSequences:
        os << "remove";
        break;
      case OccurrenceFilter::MaskOccurrences:
        os << "mask";
        break;
    }
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

