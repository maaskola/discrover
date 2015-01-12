/*
 * =====================================================================================
 *
 *       Filename:  score.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <cmath>
#include <iostream>
#include "../stats_config.hpp"
#include "correction.hpp"
#include "motif.hpp"
#include "score.hpp"

const bool excessive_debug = false;

using std::cout;
using std::endl;
using std::string;
using std::vector;

vector_t row_sums(const matrix_t &m) {
  size_t n = m.size1();
  vector_t v = boost::numeric::ublas::scalar_vector<fp_t>(n, 0);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < m.size2(); j++)
      v(i) += m(i, j);
  return v;
}

vector_t col_sums(const matrix_t &m) {
  size_t n = m.size2();
  vector_t v = boost::numeric::ublas::scalar_vector<fp_t>(n, 0);
  for (size_t j = 0; j < m.size1(); j++)
    for (size_t i = 0; i < n; i++)
      v(i) += m(j, i);
  return v;
}

double compute_mutual_information_variance(const matrix_t &m_,
                                           double pseudo_count,
                                           bool normalize) {
  matrix_t m = m_;
  double z = 1;
  if (pseudo_count != 0)
    normalize = true;
  if (normalize) {
    z = 0;
    for (size_t i = 0; i < m.size1(); i++)
      for (size_t j = 0; j < m.size2(); j++)
        z += (m(i, j) += pseudo_count);
    for (size_t i = 0; i < m.size1(); i++)
      for (size_t j = 0; j < m.size2(); j++)
        m(i, j) /= z;
  }

  double mi = 0;

  // Compute row and column sums
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);

  // Log transform the row and column sums
  std::transform(begin(rs), end(rs), begin(rs), log);
  std::transform(begin(cs), end(cs), begin(cs), log);

  // Compute mutual information
  for (size_t i = 0; i < m.size1(); i++)
    for (size_t j = 0; j < m.size2(); j++)
      if (m(i, j) != 0) {
        double x = log(m(i, j)) - rs(i) - cs(j);
        mi += m(i, j) * x * x;
      }

  // Scale to a binary logarithmic base
  mi /= log(2.0);

  double j = compute_mutual_information(m_, pseudo_count, normalize, true);
  double var = 1 / (z + 1) * (mi - j * j);
  return var;
}

double compute_mutual_information(const matrix_t &m_, double pseudo_count,
                                  bool normalize, bool do_correction) {
  matrix_t m = m_;
  double z = 1;
  if (pseudo_count != 0)
    normalize = true;
  if (normalize) {
    z = 0;
    for (size_t i = 0; i < m.size1(); i++)
      for (size_t j = 0; j < m.size2(); j++)
        z += (m(i, j) += pseudo_count);
    for (size_t i = 0; i < m.size1(); i++)
      for (size_t j = 0; j < m.size2(); j++)
        m(i, j) /= z;
  }

  double mi = 0;

  // Compute row and column sums
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);

  // Log transform the row and column sums
  std::transform(begin(rs), end(rs), begin(rs), log);
  std::transform(begin(cs), end(cs), begin(cs), log);

  // Compute mutual information
  for (size_t i = 0; i < m.size1(); i++)
    for (size_t j = 0; j < m.size2(); j++)
      if (m(i, j) != 0)
        mi += m(i, j) * (log(m(i, j)) - rs(i) - cs(j));

  // Scale to a binary logarithmic base
  mi /= log(2.0);

  if (do_correction) {
    mi += (m.size1() - 1) * (m.size2() - 1) / (z + 1);
  }
  return mi;
}

double compute_gtest(const matrix_t &m, double pseudo_count, bool normalize) {
  double n = 0;
  for (size_t i = 0; i < m.size1(); i++)
    for (size_t j = 0; j < m.size2(); j++)
      n += m(i, j);
  n += m.size1() * m.size2() * pseudo_count;
  return log(2.0) * 2 * n
         * compute_mutual_information(m, pseudo_count, normalize);
}

double compute_logp_gtest(const matrix_t &m, double pseudo_count,
                          bool normalize) {
  return -pchisq(compute_gtest(m, pseudo_count, normalize),
                 (m.size1() - 1) * (m.size2() - 1), false, true);
}

double compute_bonferroni_corrected_logp_gtest(const matrix_t &m,
                                               double pseudo_count,
                                               bool normalize, size_t length,
                                               size_t degeneracy,
                                               bool dynamic_mode) {
  double log_correction;
  if (dynamic_mode)
    log_correction = compute_correction(length, degeneracy);
  else
    log_correction = log(149);

  return compute_logp_gtest(m, pseudo_count, normalize) - log_correction;
}

double compute_mcc(double a, double b, double c, double d) {
  return (a * d - b * c) / sqrt((a + b) * (a + c) * (b + d) * (c + d));
}

double compute_delta_frequency(double a, double b, double c, double d) {
  return a / (a + b) - c / (c + d);
}

count_vector_t reduce_count(const count_vector_t &full,
                            const Seeding::Collection &collection,
                            const string &contrast_name) {
  count_vector_t counts(full.size());
  auto iter = begin(counts);
  auto full_iter = begin(full);
  for (auto &contrast : collection) {
    size_t n = contrast.sets.size();
    if (contrast.name == contrast_name) {
      std::copy(full_iter, full_iter + n, iter);
      iter += n;
      full_iter += n;
    } else
      full_iter += n;
  }
  counts.resize(std::distance(begin(counts), iter));
  return counts;
}

double compute_score(const Seeding::Collection &collection,
                     const Seeding::Result &result,
                     const Seeding::Options &options,
                     Measures::Discrete::Measure measure, bool do_correction) {
  return compute_score(
      collection, result.counts, options, result, result.motif.length(),
      Seeding::motif_degeneracy(result.motif), measure, do_correction);
}

double compute_score(const Seeding::Collection &collection,
                     const count_vector_t &counts,
                     const Seeding::Options &options,
                     const Seeding::Objective &objective, size_t length,
                     size_t degeneracy, Measures::Discrete::Measure measure,
                     bool do_correction) {
  // double compute_score(const Seeding::Collection &collection, const
  // count_vector_t &counts, const Seeding::Options &options,
  // Measures::Discrete::Measure measure, size_t length, size_t degeneracy, bool
  // do_correction) {
  if (options.verbosity >= Verbosity::debug) {
    cout << "compute_score(Seeding::Collection)" << endl;
    cout << "counts = " << vec2string(counts) << endl;
  }
  if (measure == Measures::Discrete::Measure::Undefined)
    measure = objective.measure;
  if (options.verbosity >= Verbosity::debug)
    cout << "Using measure = " << measure2string(measure) << endl;
  double score = 0;
  double W = 0;
  for (auto &expr : objective) {
    auto contrast_iter = collection.find(expr.contrast);
    if (contrast_iter == end(collection)) {
      vector<string> names;
      for (auto &x : collection)
        names.push_back(x.name);
      throw Exception::Plasma::NoContrastForObjective(to_string(expr), names);
    } else {
      auto counts_ = reduce_count(counts, collection, expr.contrast);
      double w = options.weighting ? contrast_iter->set_size : 1;
      W += w;
      if (options.verbosity >= Verbosity::debug)
        cout << "counts_ = " << vec2string(counts_) << endl;
      score
          += expr.sign * w
             * compute_score(*contrast_iter, counts_, options, measure, length,
                             degeneracy, objective.motif_name, do_correction);
    }
  }
  if (options.weighting)
    score /= W;
  if (options.verbosity >= Verbosity::debug)
    cout << "end: compute_score(Seeding::Collection) -> " << score << endl;
  return score;
}

double compute_score(const Seeding::Contrast &contrast,
                     const count_vector_t &counts,
                     const Seeding::Options &options,
                     Measures::Discrete::Measure measure, size_t length,
                     size_t degeneracy, const string &motif_name,
                     bool do_correction) {
  if (measure == Measures::Discrete::Measure::Undefined)
    throw Exception::Plasma::UndefinedMeasure();
  if (measure != Measures::Discrete::Measure::SignalFrequency
      and contrast.sets.size() < 2)
    return -std::numeric_limits<double>::infinity();

  matrix_t occurrence_table = matrix_t(counts.size(), 2);
  for (size_t i = 0; i < counts.size(); i++) {
    occurrence_table(i, 0) = counts[i];
    occurrence_table(i, 1) = (options.word_stats ? contrast.sets[i].seq_size
                                                 : contrast.sets[i].set_size)
                             - counts[i];
  }
  double signal_freq = 0, control_freq = 0;
  size_t signal_size = 0, control_size = 0;

  size_t sample_idx = 0;
  bool signal_present = false;
  for (auto &dataset : contrast) {
    if (excessive_debug and options.verbosity >= Verbosity::debug) {
      std::cerr << "motifs in dataset:";
      for (auto &x : dataset.motifs)
        std::cerr << " " << x;
      std::cerr << std::endl;
      std::cerr << "looking for: " << motif_name;
    }
    if (find(begin(dataset.motifs), end(dataset.motifs), motif_name)
        != end(dataset.motifs)) {
      signal_freq += occurrence_table(sample_idx, 0);
      signal_size += (options.word_stats ? dataset.seq_size : dataset.set_size);
      signal_present = true;
      if (excessive_debug and options.verbosity >= Verbosity::debug)
        std::cerr << " - found!" << std::endl;
    } else {
      control_freq += occurrence_table(sample_idx, 0);
      control_size
          += (options.word_stats ? dataset.seq_size : dataset.set_size);
      if (excessive_debug and options.verbosity >= Verbosity::debug)
        std::cerr << " - not found!" << std::endl;
    }
    sample_idx++;
  }
  double signal_rel_freq = signal_freq / signal_size;
  double control_rel_freq = control_freq / control_size;

  if (excessive_debug and options.verbosity >= Verbosity::debug)
    std::cerr << "options.no_enrichment_filter = "
              << options.no_enrichment_filter
              << " signal_present = " << signal_present
              << " signal_rel_freq = " << signal_rel_freq
              << " control_rel_freq = " << control_rel_freq << std::endl;
  if (not options.no_enrichment_filter and signal_present
      and signal_rel_freq < control_rel_freq)
    return -std::numeric_limits<double>::infinity();

  double score = 0;
  switch (measure) {
    case Measures::Discrete::Measure::MutualInformation:
      score = compute_mutual_information(occurrence_table, options.pseudo_count,
                                         true, do_correction);
      break;
    case Measures::Discrete::Measure::VarianceMutualInformation:
      score = compute_mutual_information_variance(occurrence_table,
                                                  options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::Gtest:
      score = compute_gtest(occurrence_table, options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::LogpGtest:
      score = compute_logp_gtest(occurrence_table, options.pseudo_count, true);
      break;
    case Measures::Discrete::Measure::CorrectedLogpGtest:
      score = compute_bonferroni_corrected_logp_gtest(
          occurrence_table, options.pseudo_count, true, length, degeneracy,
          not options.fixed_motif_space_mode);
      break;
    case Measures::Discrete::Measure::MatthewsCorrelationCoefficient:
      return compute_mcc(signal_freq, signal_size - signal_freq, control_freq,
                         control_size - control_freq);
    case Measures::Discrete::Measure::DeltaFrequency:
      return signal_rel_freq - control_rel_freq;
    case Measures::Discrete::Measure::SignalFrequency:
      return signal_rel_freq;
    case Measures::Discrete::Measure::ControlFrequency:
      return control_rel_freq;
    default:
      assert(false);
  }

  return score;
}

namespace Exception {
namespace Plasma {
const char *UndefinedMeasure::what() const noexcept {
  string msg
      = "Error: undefined score (Measures::Discrete::Measure::undefined).";
  return msg.c_str();
}
NoContrastForObjective::NoContrastForObjective(
    const std::string &expr_, const std::vector<std::string> &names_)
    : exception(), expr(expr_), names(names_){};
const char *NoContrastForObjective::what() const noexcept {
  string msg
      = "Error: couldn't find contrast for objective contrast expression: '"
        + expr + "'.\n";
  for (auto &name : names)
    msg += "\nThis contrast is present: " + name;
  return msg.c_str();
}
}
}
