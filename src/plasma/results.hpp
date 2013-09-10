/*
 * =====================================================================================
 *
 *       Filename:  results.hpp
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


#ifndef RESULTS_HPP
#define RESULTS_HPP

#include <string>
#include "stats.hpp"
#include "options.hpp"

namespace Plasma {
  struct Result : public Objective {
    std::string motif;
    double score;
    double log_p;
    Stats::OccurrenceCounts counts;
    Result(const Objective &objective);
  };
  typedef std::vector<Result> Results;
}

#endif   /* ----- #ifndef RESULTS_HPP  ----- */

