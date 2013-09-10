/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  A tool MIC
 *
 *        Version:  1.0
 *        Created:  08/03/2012 04:30:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "association.hpp"
#include "verbosity.hpp"

namespace MIC {

  typedef double score_t;
  //typedef size_t data_t;
  typedef double data_t;
  typedef std::vector<data_t> sample_t;
  typedef std::vector<size_t> bounds_t;
  typedef sample_t stats_t;

  bounds_t sample(const sample_t &d, const std::set<size_t> &fixed_bounds, size_t n, size_t n_iter, Verbosity verbosity);

};

