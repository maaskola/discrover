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
#include <set>
#include <vector>
#include "association.hpp"
#include "verbosity.hpp"
#include "mic_impl.hpp"

namespace MIC {

  std::vector<data_t> read(const std::string &path)
  {
    std::vector<data_t> v;
    std::ifstream ifs(path.c_str());
    std::string line;
    while(getline(ifs, line)) {
      data_t x = atof(line.c_str());
      v.push_back(x);
    }
    return(v);
  }

  void doit_simple(const std::vector<size_t> &v)
  {
    for(size_t i = 1; i < v.size(); i++) {
      matrix_t m = zero_matrix(2,2);
      for(size_t j = 0; j < i; j++)
        if(v[j] == 1)
          m(0,0)++;
        else
          m(0,1)++;
      for(size_t j = i; j < v.size(); j++)
        if(v[j] == 1)
          m(1,0)++;
        else
          m(1,1)++;
      double mi = calc_mutual_information(m, 0, true, false, false);
      std::cout << "i = " << i << " mi = " << mi << " mic = " << mi / log(2) << std::endl;
    }
  }


  void get_stats(const sample_t &d, const bounds_t &ends, stats_t &stats)
  {
    size_t last = 0;
    for(size_t i = 0; i < ends.size(); i++) {
      for(size_t j = last; j < ends[i]; j++)
        stats[i] += d[j];
      last = ends[i];
    }
  }

  void print_state(const bounds_t &ends, stats_t &stats)
  {
    size_t last = 0;
    for(size_t i = 0; i < ends.size(); i++) {
      std::cout << i << "\t" << ends[i] - last << ":" << last << ".." << ends[i] << "\t" << stats[i] << "\t" << stats[i] / (ends[i] - last) << std::endl;
      last = ends[i];
    }
  }

  score_t calc_score(const bounds_t &ends, const stats_t &stats, Verbosity verbosity)
  {
    size_t k = stats.size();
    matrix_t m = zero_matrix(k,2);
    size_t last = 0;
    for(size_t j = 0; j < k; j++) {
      m(j,0) = stats[j];
      m(j,1) = (ends[j] - last) - stats[j];
      last = ends[j];
    }
    double mi = calc_mutual_information(m, 1, true, false, false);
    if(verbosity >= Verbosity::verbose)
      std::cout << "mi = " << mi << std::endl;
    // std::cout << "i = " << i << " mi = " << mi << " mic = " << mi / log(2) << std::endl;
    return(mi);
  }


  bool acceptable(score_t cur, score_t prev, double T, Verbosity verbosity)
  {
    bool ok = cur >= prev;
    if(verbosity >= Verbosity::verbose)
      std::cout << "Better? " << (ok ? "yes!" : "no!") << std::endl;
    if(not ok) {
      double d = cur - prev;
      double r = 1.0 * rand() / RAND_MAX;
      double p = exp(d*T);
      if(r <= p)
        ok = true;
      if(verbosity >= Verbosity::verbose)
        std::cout << "Acceptable? " << (ok ? "yes!" : "no!") << " r = " << r << " p = " << p << " d = " << d << " T = " << T << std::endl;
    }
    return(ok);
  }

  void try_modification(const std::set<size_t> &fixed, const sample_t &d, size_t n, size_t l, bounds_t &ends, stats_t &stats, double &score, double &temperature, Verbosity verbosity)
  {
    bounds_t new_ends(ends);
    stats_t new_stats(stats);

    // select one end to remove
    size_t selected = rand() % (n-1);
    while(fixed.find(ends[selected]) != fixed.end()) // do not choose one of the fixed boundaries
      selected = rand() % (n-1);

    // determine a new end position
    size_t new_end = 1 + (rand() % (l-1));

    // choose a new position until there's not already an end at this position
    while(std::find(ends.begin(), ends.end(), new_end) != ends.end())
      new_end = 1 + (rand() % (l-1));

    if(verbosity >= Verbosity::verbose) {
      std::cout << "selected the " << selected << "th end to be removed" << std::endl;
      std::cout << "chose " << new_end << " as a new end position" << std::endl;
    }

    // move stats of selected end to the following end
    new_stats[selected+1] += new_stats[selected];

    // remove the selected end and its associated stats
    new_stats.erase(new_stats.begin() + selected);
    new_ends.erase(new_ends.begin() + selected);

    // determine which stats have to be updated
    selected = 0;
    while(new_end >= new_ends[selected])
      selected++;

    data_t new_count = 0;
    size_t first = 0;
    if(selected > 0)
      first = new_ends[selected-1];
    for(size_t i = first; i < new_end; i++)
      new_count += d[i];

    new_stats[selected] -= new_count;
    new_stats.insert(new_stats.begin() + selected, new_count);
    new_ends.insert(new_ends.begin() + selected, new_end);

    if(verbosity >= Verbosity::verbose) {
      std::cout << "new_count = " << new_count << std::endl;
      std::cout << "Candidate state:" << std::endl;
      print_state(new_ends, new_stats);
    }

    score_t new_score = calc_score(new_ends, new_stats, verbosity);
    if(acceptable(new_score, score, temperature, verbosity)) {
      if(verbosity >= Verbosity::verbose)
        std::cout << "Accepting new state" << std::endl;
      score = new_score;
      stats = new_stats;
      ends = new_ends;
    }
  }

  bounds_t sample(const sample_t &d, const std::set<size_t> &fixed_bounds, size_t n_, size_t n_iter, Verbosity verbosity)
  {
    size_t n = n_ + fixed_bounds.size();
    size_t l = d.size();
    if(verbosity >= Verbosity::verbose)
      std::cout << "Total length = " << l << std::endl;
    bounds_t ends;
    for(auto &b: fixed_bounds)
      ends.push_back(b);
    if(l + 1 < n)
      return(ends);
    if(l == n) {
      bounds_t ends(n);
      for(size_t i = 0; i < n; i++)
        ends[i] = i+1;
      return(ends);
    }
    for(size_t i = 0; i < n_; i++) {
      size_t k = rand() % l;
      while(std::find(ends.begin(), ends.end(), k) != ends.end())
        k = rand() % l;
      ends.push_back(k);
    }
    std::sort(ends.begin(), ends.end());
    stats_t counts(n,0);
    get_stats(d, ends, counts);
    if(verbosity >= Verbosity::verbose) {
      std::cout << "Initial:" << std::endl;
      print_state(ends, counts);
    }

    score_t score = calc_score(ends, counts, verbosity);

    bounds_t best_ends(ends);
    stats_t best_counts(counts);
    double best_score = score;

    double temperature = 1;
    while(n_iter-- > 0) {
      try_modification(fixed_bounds, d, n, l, ends, counts, score, temperature, verbosity);
      if(score > best_score) {
        best_score = score;
        best_ends = ends;
        best_counts = counts;
      }
      temperature += 0.02;
      if(verbosity >= Verbosity::verbose) {
        std::cout << "Current state:" << std::endl;
        print_state(ends, counts);
      }
    }

    if(verbosity >= Verbosity::info) {
      std::cout << "Final:" << std::endl;
      print_state(best_ends, best_counts);
      std::cout << "MIC: " << n << " " << calc_score(best_ends, best_counts, verbosity) << std::endl;
    }
    return(best_ends);
  }


  /*
     using namespace std;
     int main(int argc, const char** argv)
     {
     srand(time(0));
     vector<data_t> v = { 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
     if(argc == 2)
     v = read(argv[1]);
  // doit_simple(v);
  // bounds_t bounds = sample(v, 20, 1000000, Verbosity::info);
  data_t n = 0;
  for(auto x: v)
  n += x;
  double limit = pow(n, 0.6);
  limit = std::min<double>(5, limit);
  cout << "Limit = " << limit << endl;
  for(size_t i = 2; i <= limit; i++)
  bounds_t bounds = sample(v, i, 1000000, Verbosity::info);
  return(EXIT_SUCCESS);
  }
  */
};

