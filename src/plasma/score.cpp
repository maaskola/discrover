/*
 * =====================================================================================
 *
 *       Filename:  score.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <cmath>
#include <iostream>
#include "../stats_config.hpp"
#include "correction.hpp"
#include "score.hpp"

const bool excessive_debug = false;

using std::cout;
using std::endl;
using std::string;

vector_t row_sums(const matrix_t &m) {
  size_t n = m.size1();
  vector_t v = boost::numeric::ublas::scalar_vector<fp_t>(n,0);
  for(size_t i = 0; i < n; i++)
    for(size_t j = 0; j < m.size2(); j++)
      v(i) += m(i,j);
  return(v);
}

vector_t col_sums(const matrix_t &m) {
  size_t n = m.size2();
  vector_t v = boost::numeric::ublas::scalar_vector<fp_t>(n,0);
  for(size_t j = 0; j < m.size1(); j++)
    for(size_t i = 0; i < n; i++)
      v(i) += m(j,i);
  return(v);
}

double compute_mutual_information_variance(const Plasma::Stats::OccurrenceTable &m_, double pseudo_count, bool normalize) {
  Plasma::Stats::OccurrenceTable m = m_;
  double z = 1;
  if(pseudo_count != 0)
    normalize = true;
  if(normalize) {
    z = 0;
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        z += (m(i,j) += pseudo_count);
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        m(i,j) /= z;
  }

  double mi = 0;

  // Compute row and column sums
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);

  // Log transform the row and column sums
  std::transform(begin(rs), end(rs), begin(rs), log);
  std::transform(begin(cs), end(cs), begin(cs), log);

  // Compute mutual information
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      if(m(i,j) != 0) {
        double x = log(m(i,j)) - rs(i) - cs(j);
        mi += m(i,j) * x * x;
      }

  // Scale to a binary logarithmic base
  mi /= log(2.0);

  double j = compute_mutual_information(m_, pseudo_count, normalize, true);
  double var = 1/ (z+1) * (mi - j*j);
  return(var);
}

double compute_mutual_information(const Plasma::Stats::OccurrenceTable &m_, double pseudo_count, bool normalize, bool do_correction) {
  Plasma::Stats::OccurrenceTable m = m_;
  double z = 1;
  if(pseudo_count != 0)
    normalize = true;
  if(normalize) {
    z = 0;
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        z += (m(i,j) += pseudo_count);
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        m(i,j) /= z;
  }

  double mi = 0;

  // Compute row and column sums
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);

  // Log transform the row and column sums
  std::transform(begin(rs), end(rs), begin(rs), log);
  std::transform(begin(cs), end(cs), begin(cs), log);

  // Compute mutual information
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      if(m(i,j) != 0)
        mi += m(i,j) * (log(m(i,j)) - rs(i) - cs(j));

  // Scale to a binary logarithmic base
  mi /= log(2.0);

  if(do_correction) {
    mi += (m.size1() - 1) * (m.size2() - 1) / (z + 1);
  }
  return(mi);
}

double compute_gtest(const Plasma::Stats::OccurrenceTable &m, double pseudo_count, bool normalize) {
  double n = 0;
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      n += m(i,j);
  n += m.size1() * m.size2() * pseudo_count;
  return(compute_mutual_information(m, pseudo_count, normalize) * log(2.0) * 2 * n);
}

double compute_logp_gtest(const Plasma::Stats::OccurrenceTable &m, double pseudo_count, bool normalize) {
  return(-pchisq(compute_gtest(m, pseudo_count, normalize), (m.size1() - 1) * (m.size2() - 1), false, true));
}

double compute_bonferroni_corrected_logp_gtest(const Plasma::Stats::OccurrenceTable &m, double pseudo_count, bool normalize, size_t length, size_t degeneracy) {
  double log_correction = compute_correction(length, degeneracy);
  return(compute_logp_gtest(m, pseudo_count, normalize) - log_correction);
}

double compute_mcc(double a, double b, double c, double d) {
  return (a*d - b*c) / sqrt((a+b)*(a+c)*(b+d)*(c+d));
}

double compute_delta_frequency(double a, double b, double c, double d) {
  return a / (a + b) - c / (c + d);
}

Plasma::Stats::OccurrenceCounts reduce_count(const Plasma::Stats::OccurrenceCounts &full, const Plasma::DataCollection &collection, const string &series_name) {
  Plasma::Stats::OccurrenceCounts counts(full.size());
  auto iter = begin(counts);
  auto full_iter = begin(full);
  for(auto &series: collection) {
    size_t n = series.sets.size();
    if(series.name == series_name) {
      std::copy(full_iter, full_iter+n, iter);
      iter += n;
      full_iter += n;
    } else
      full_iter += n;
  }
  counts.resize(std::distance(begin(counts), iter));
  return(counts);
}

double compute_score(const Plasma::DataCollection &collection, const Plasma::Result &result, const Plasma::options_t &options, Measures::Discrete::Measure measure, bool do_correction) {
  return(compute_score(collection, result.counts, options, result, result.motif.length(), Plasma::motif_degeneracy(result.motif), measure, do_correction));
}

double compute_score(const Plasma::DataCollection &collection,
    const Plasma::Stats::OccurrenceCounts &counts,
    const Plasma::options_t &options,
    const Plasma::Objective &objective,
    size_t length,
    size_t degeneracy,
    Measures::Discrete::Measure measure,
    bool do_correction) {
// double compute_score(const Plasma::DataCollection &collection, const Plasma::Stats::OccurrenceCounts &counts, const Plasma::options_t &options, Measures::Discrete::Measure measure, size_t length, size_t degeneracy, bool do_correction) {
  if(options.verbosity >= Verbosity::debug) {
    cout << "compute_score(Plasma::DataCollection)" << endl;
    cout << "counts = " << counts << endl;
  }
  if(measure == Measures::Discrete::Measure::Undefined)
    measure = objective.measure;
  if(options.verbosity >= Verbosity::debug)
    cout << "Using measure = " << measure2string(measure) << endl;
  double score = 0;
  for(auto &expr: objective) {
    auto series_iter = collection.find(expr.series);
    if(series_iter == end(collection)) {
      cout << "Error: couldn't find series for objective series expression: " << expr << endl;
      for(auto &x: collection)
        cout << "This series is present: " << x.name << std::endl;
      exit(-1);
    } else {
      auto counts_ = reduce_count(counts, collection, expr.series);
      if(options.verbosity >= Verbosity::debug)
        cout << "counts_ = " << counts_ << endl;
      score += expr.sign * compute_score(*series_iter, counts_, options, measure, length, degeneracy, objective.motif_name, do_correction);
    }
  }
  if(options.verbosity >= Verbosity::debug)
    cout << "end: compute_score(Plasma::DataCollection) -> " << score << endl;
  return(score);
}

double compute_score(const Plasma::DataSeries &data_series, const Plasma::Stats::OccurrenceCounts &counts, const Plasma::options_t &options, Measures::Discrete::Measure measure, size_t length, size_t degeneracy, const string &motif_name, bool do_correction) {
  if(measure == Measures::Discrete::Measure::Undefined) {
    cout << "Error in compute_score: Measures::Discrete::Measure::undefined." << endl;
    exit(-1);
  }
  if(measure != Measures::Discrete::Measure::SignalFrequency and data_series.sets.size() < 2)
    return(-std::numeric_limits<double>::infinity());

  matrix_t occurrence_table = matrix_t(counts.size(), 2);
  for(size_t i = 0; i < counts.size(); i++) {
    occurrence_table(i,0) = counts[i];
    occurrence_table(i,1) = (options.word_stats ? data_series.sets[i].seq_size : data_series.sets[i].set_size) - counts[i];
  }
  double signal_freq = 0, control_freq = 0;
  size_t signal_size = 0, control_size = 0;

  size_t sample_idx = 0;
  bool signal_present = false;
  for(auto &data_set: data_series) {
    if(excessive_debug and options.verbosity >= Verbosity::debug) {
      std::cerr << "motifs in data_set:";
      for(auto &x: data_set.motifs)
        std::cerr << " " << x;
      std::cerr << std::endl;
      std::cerr << "looking for: " << motif_name;
    }
    if(find(begin(data_set.motifs), end(data_set.motifs), motif_name) != end(data_set.motifs)) {
      signal_freq += occurrence_table(sample_idx, 0);
      signal_size += (options.word_stats ? data_set.seq_size : data_set.set_size);
      signal_present = true;
      if(excessive_debug and options.verbosity >= Verbosity::debug)
        std::cerr << " - found!" << std::endl;
    } else {
      control_freq += occurrence_table(sample_idx, 0);
      control_size += (options.word_stats ? data_set.seq_size : data_set.set_size);
      if(excessive_debug and options.verbosity >= Verbosity::debug)
        std::cerr << " - not found!" << std::endl;
    }
    sample_idx++;
  }
  double signal_rel_freq = signal_freq / signal_size;
  double control_rel_freq = control_freq / control_size;

  if(excessive_debug and options.verbosity >= Verbosity::debug)
    std::cerr << "options.no_enrichment_filter = " << options.no_enrichment_filter
      << " signal_present = " << signal_present
      << " signal_rel_freq = " << signal_rel_freq
      << " control_rel_freq = " << control_rel_freq << std::endl;
  if(not options.no_enrichment_filter and signal_present and signal_rel_freq < control_rel_freq)
    return(-std::numeric_limits<double>::infinity());

  double score = 0;
  switch(measure) {
    case Measures::Discrete::Measure::MutualInformation:
      score = compute_mutual_information(occurrence_table, options.pseudo_count, true, do_correction);
      break;
    case Measures::Discrete::Measure::VarianceMutualInformation:
      score = compute_mutual_information_variance(occurrence_table, options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::Gtest:
      score = compute_gtest(occurrence_table, options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::LogpGtest:
      score = compute_logp_gtest(occurrence_table, options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::CorrectedLogpGtest:
      score = compute_bonferroni_corrected_logp_gtest(occurrence_table, options.pseudo_count, true, length, degeneracy);
      break;
    case Measures::Discrete::Measure::MatthewsCorrelationCoefficient:
      return(compute_mcc(signal_freq, signal_size - signal_freq, control_freq, control_size - control_freq));
    case Measures::Discrete::Measure::DeltaFrequency:
      return(signal_rel_freq - control_rel_freq);
    case Measures::Discrete::Measure::SignalFrequency:
      return(signal_rel_freq);
    case Measures::Discrete::Measure::ControlFrequency:
      return(control_rel_freq);
    default:
      assert(false);
  }

  return(score);
}
