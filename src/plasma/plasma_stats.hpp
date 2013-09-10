/*
 * =====================================================================================
 *
 *       Filename:  plasma_stats.cpp
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


#ifndef  PLASMA_STATS_HPP
#define  PLASMA_STATS_HPP

#include <unordered_map>
#include <map>
#include <string>
#include "plasma.hpp"
#include "motif.hpp"
#include "stats.hpp"

namespace Plasma {
  typedef std::unordered_map<std::string, Stats::OccurrenceCounts> hash_map_t;
  typedef std::unordered_map<std::string, double> score_map_t;
  typedef std::multimap<double, std::string> rev_map_t;
}

#endif   /* ----- #ifndef PLASMA_STATS_HPP ----- */

