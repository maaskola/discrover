/* =====================================================================================
 * Copyright (c) 2011, Jonas Maaskola
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 *
 *       Filename:  hmm_gradient.cpp
 *
 *    Description:  Gradient calculation routines of the HMM class 
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <omp.h>
#include <cmath>
#include "../aux.hpp"
#include "hmm.hpp"
#include "logistic.hpp"
#include "subhmm.hpp"

using namespace std;

#define DO_PARALLEL 1

double HMM::log_likelihood_gradient(const Data::Seqs &seqs, const Training::Targets &targets, matrix_t &transition_g, matrix_t &emission_g) const
{
  transition_g = zero_matrix(n_states, n_states);
  emission_g = zero_matrix(n_states, n_emissions);

  double lp = 0;
  vector<matrix_t> t_g, e_g;

#pragma omp parallel shared(emission_g, transition_g) if(DO_PARALLEL)
  {
    // Initalize storage for thread intermediate results
#pragma omp single
    {
      size_t n_threads = omp_get_num_threads();
      //    cout << "B Num threads = " << n_threads << endl;
      if(not targets.transition.empty())
        t_g = vector<matrix_t> (n_threads, zero_matrix(transition_g.size1(), transition_g.size2()));
      if(not targets.emission.empty())
        e_g = vector<matrix_t> (n_threads, zero_matrix(emission_g.size1(), emission_g.size2()));
    }

    // Compute likelihood for each sequence
#pragma omp for reduction (+:lp)
    for(size_t i = 0; i < seqs.size(); i++) {
      int thread_idx = omp_get_thread_num();
      vector_t scale;
      matrix_t f = compute_forward_scaled(seqs[i], scale);
      matrix_t b = compute_backward_prescaled(seqs[i], scale);

      // Compute expected statistics
      matrix_t T, E;
      double logp = BaumWelchIteration_single(T, E, seqs[i], targets);

      if(not targets.transition.empty())
        // Compute log likelihood gradients w.r.t. transition probability 
        t_g[thread_idx] += transition_gradient(T, targets.transition);

      if(not targets.emission.empty())
        // Compute log likelihood gradients w.r.t. emission probability
        e_g[thread_idx] += emission_gradient(E, targets.emission);

      lp += logp;
    }

    // Collect results of threads
#pragma omp single
    {
      if(not targets.transition.empty())
        for(auto &x: t_g)
          transition_g += x;
      if(not targets.emission.empty())
        for(auto &x: e_g)
          emission_g += x;
    }
  }

  return(lp);
}

double HMM::chi_square_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::chi_square_gradient(Data::Contrast, Feature)" << endl;
  assert(0);
  size_t n_samples = contrast.sets.size();
  matrix_t counts(n_samples, 2);
  vector<matrix_t> trans_g(n_samples), emission_g(n_samples);

  if(verbosity >= Verbosity::debug)
    cout << "Looking at groups " << present << endl;

  for(size_t sample_idx = 0; sample_idx < n_samples; sample_idx++) {
    counts(sample_idx, 0) = posterior_gradient(contrast.sets[sample_idx], task, present, trans_g[sample_idx], emission_g[sample_idx]).posterior;
    counts(sample_idx, 1) = contrast.sets[sample_idx].set_size - counts(sample_idx, 0);
  }

  counts = counts + pseudo_count; // Add pseudo-count

  if(verbosity >= Verbosity::debug)
    cout << "comparison: " << counts << endl;

  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  double total_r = 0;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++) {
    total_r += counts(set_idx, 1);
    double log_r = log(counts(set_idx, 1));

    if(verbosity >= Verbosity::debug)
      cout << "log_r = " << log_r << endl;

    double factor = log(counts(set_idx, 0)) - log_r;

    if(verbosity >= Verbosity::debug)
      cout << "set_idx = " << set_idx << " factor = " << factor << endl;
    if(not task.targets.transition.empty())
      g.transition += factor * trans_g[set_idx];
    if(not task.targets.emission.empty())
      g.emission += factor * emission_g[set_idx];
  }
  double log_total_r = log(total_r);

  double sum_post = 0;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++)
    sum_post += counts(set_idx, 0);

  double factor = log(sum_post) - log_total_r;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++) {
    if(not task.targets.transition.empty())
      g.transition -= factor * trans_g[set_idx];
    if(not task.targets.emission.empty())
      g.emission -= factor * emission_g[set_idx];
  }

  const bool divide_by_data_size = false;
  if(divide_by_data_size) {
    double n = contrast.set_size + pseudo_count * counts.size1() * counts.size2();
    if(not task.targets.transition.empty())
      g.transition /= n;
    if(not task.targets.emission.empty())
      g.emission /= n;
  }

  if(verbosity >= Verbosity::verbose) {
    if(not task.targets.transition.empty())
      cout << "Summary posterior transition gradient = " << g.transition << endl;
    if(not task.targets.emission.empty())
      cout << "Summary posterior emission gradient = " << g.emission << endl;
  }

  double mi = calc_mutual_information(counts, 0, true, false, false);
  if(verbosity >= Verbosity::debug)
    cout << "evaluation: " << mi << endl;
  return(mi);
}

/** Gradient of Matthew's correlation coefficient. */
double HMM::matthews_correlation_coefficient_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information_gradient(Data::Contrast, Feature)" << endl;
  size_t n_samples = contrast.sets.size();

  Gradient signal_g, control_g;
  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.transition.empty()) {
    signal_g.transition = zero_matrix(n_states, n_states);
    control_g.transition = zero_matrix(n_states, n_states);
  }
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);
  if(not task.targets.emission.empty()) {
    signal_g.emission = zero_matrix(n_states, n_emissions);
    control_g.emission = zero_matrix(n_states, n_emissions);
  }

  double signal_counts = 0, control_counts = 0;
  size_t n_signal = 0, n_control = 0;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < n_samples; set_idx++) {
    matrix_t trans_g = zero_matrix(n_states, n_states);
    matrix_t emission_g = zero_matrix(n_states, n_emissions);

    double counts = posterior_gradient(contrast.sets[set_idx], task, present, trans_g, emission_g).posterior;

    bool signal = is_present(contrast.sets[set_idx], present);
    size_t total = contrast.sets[set_idx].set_size;
    total += 2 * pseudo_count; // Add pseudo-count
    counts += pseudo_count; // Add pseudo-count
    if(signal) {
      n_signal += total;
      signal_counts += counts;
      if(not task.targets.transition.empty())
        signal_g.transition += trans_g;
      if(not task.targets.emission.empty())
        signal_g.emission += emission_g;
    } else {
      n_control += total;
      control_counts += counts;
      if(not task.targets.transition.empty())
        control_g.transition += trans_g;
      if(not task.targets.emission.empty())
        control_g.emission += emission_g;
    }

    if(verbosity >= Verbosity::debug)
      cout << "set_idx = " << set_idx << " path = " << contrast.sets[set_idx].path << " signal = " << signal << " counts = " << counts << endl
        << "current emission gradient = " << emission_g << endl;
  }

  confusion_matrix counts;
  counts.true_positives = signal_counts;
  counts.false_positives = control_counts;
  counts.false_negatives = n_signal - signal_counts;
  counts.true_negatives = n_control - control_counts;

  if(verbosity >= Verbosity::debug) {
    cout << "confusion matrix: TP " << counts.true_positives << " FN " << counts.false_negatives << " FP " << counts.false_positives << " TN " << counts.true_negatives << endl;
    cout << "signal emission gradient = " << signal_g.emission << endl;
    cout << "control emission gradient = " << control_g.emission << endl;
  }

  double total = n_signal + n_control;
  double common_factor = 1.0 / sqrt(n_signal * n_control * (signal_counts + control_counts) * (total - signal_counts - control_counts));
  double marginal_factor = (n_control * signal_counts - n_signal * control_counts) * (total / 2.0 - signal_counts - control_counts) / (signal_counts + control_counts) / (total - signal_counts - control_counts);

  if(not task.targets.transition.empty())
    g.transition = common_factor * (
        n_control * signal_g.transition - n_signal * control_g.transition
        - marginal_factor * (signal_g.transition + control_g.transition));
  if(not task.targets.emission.empty())
    g.emission = common_factor * (
        n_control * signal_g.emission - n_signal * control_g.emission
        - marginal_factor * (signal_g.emission + control_g.emission));

  if(verbosity >= Verbosity::verbose) {
    if(not task.targets.transition.empty())
      cout << "Summary MCC transition gradient = " << g.transition << endl;
    if(not task.targets.emission.empty())
      cout << "Summary MCC emission gradient = " << g.emission << endl;
  }

  double mcc = calc_matthews_correlation_coefficient(counts);
  return(mcc);
}

double HMM::log_likelihood_difference_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::log_likelihood_difference_gradient(Data::Contrast, Feature)" << endl;
  size_t n_samples = contrast.sets.size();

  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  double log_likelihood_difference = 0;
  // for each of the samples
  for(size_t set_idx = 0; set_idx < n_samples; set_idx++) {
    matrix_t trans_g = zero_matrix(n_states, n_states);
    matrix_t emission_g = zero_matrix(n_states, n_emissions);
    double log_likel = log_likelihood_gradient(contrast.sets[set_idx].sequences, task.targets, trans_g, emission_g);

    double sign = is_present(contrast.sets[set_idx], present) ? 1.0 : -1.0;
    log_likelihood_difference += sign * log_likel;

    if(verbosity >= Verbosity::debug)
      cout << "set_idx = " << set_idx << " sign = " << sign << " likelihood = " << log_likel << endl;

    if(not task.targets.transition.empty())
      g.transition += sign * trans_g;
    if(not task.targets.emission.empty())
      g.emission += sign * emission_g;
  }
  return(log_likelihood_difference);
}

double HMM::class_likelihood_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::verbose)
    cout << "HMM::class_likelihood_gradient(Data::Contrast, Feature)" << endl;

  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  double l = 0;
  for(auto &dataset: contrast)
    l += class_likelihood_gradient(dataset, task, present, g);

  if(verbosity >= Verbosity::verbose)
    cout << "Data::Contrast l = " << l << endl;

  /*
  if(false) {
    double l2 = class_likelihood(contrast, group_idx);
    if(fabs(l - l2) > 1e-6) {
      cout << "Error: difference between class likelihood calculations: l " << l << " l2 = " << l2 << endl;
      exit(-1);
    } else
      cout << "OK!" << endl;
  } */
  return(l);
}

/** Compute the probability of correct classification, also known as class
 * likelihood, or maximum mutual information estimation (MMIE). This routine
 * actually only computes the transition and emission probability gradient, the
 * class parameters are reestimated with another routine. */
double HMM::class_likelihood_gradient(const Data::Set &dataset, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::verbose)
    cout << "HMM::class_likelihood_gradient(Data::Set, Feature)" << endl;
  // TODO reactivate
  cout << "Error: class_likelihood_gradient is currently not implemented." << endl;
  exit(-1);
  /*
  const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
  const double class_cond = get_class_motif_prior(dataset.sha1, group_idx);
  const double current_class_prior = get_class_prior(dataset.sha1);
  const double log_class_prior = log(current_class_prior);

  if(verbosity >= Verbosity::verbose)
    cout << "marginal_motif_prior = " << marginal_motif_prior << endl
      << "class_cond = " << class_cond << endl
      << "current_class_prior = " << current_class_prior << endl
      << "log_class_prior = " << log_class_prior << endl;

  double l = 0;
  vector<matrix_t> t_g, e_g;
#pragma omp parallel shared(g, t_g, e_g) if(DO_PARALLEL)
  {
    // Initialize storage for thread intermediate results
#pragma omp single
    {
      size_t n_threads = omp_get_num_threads();
      if(not task.targets.transition.empty())
        t_g = vector<matrix_t> (n_threads, zero_matrix(g.transition.size1(), g.transition.size2()));
      if(not task.targets.emission.empty())
        e_g = vector<matrix_t> (n_threads, zero_matrix(g.emission.size1(), g.emission.size2()));
    }

#pragma omp for schedule(static) reduction(+:l)
    for(size_t i = 0; i < dataset.set_size; i++) {
      int thread_idx = omp_get_thread_num();
      */

      /* c                     Class 1
       * not C                 Class 2
       * current_class_prior   P(c)
       * log_class_prior       log P(c)
       * p                     P(m|X)
       * class_cond            P(m|c)
       * marginal_motif_prior  P(m)
       * x                     log P(C|X)
       * exp(-x)               1 / P(C|X)
       */

      /*
      matrix_t t, e;
      posterior_gradient_t res = posterior_gradient(dataset.sequences[i], task, group_idx, t, e);
      double p = res.posterior;
      double x = 0;
      if(log_class_prior != 0)
        x = log_class_prior + log(p * class_cond / marginal_motif_prior +
            (1-p) * (1-class_cond) / (1-marginal_motif_prior));
          // (marginal_motif_prior < 1 ?  (1-p) * (1-class_cond) / (1-marginal_motif_prior) : 0));
      if(verbosity >= Verbosity::verbose)
        cout << "Sequence " << dataset.sequences[i].definition << " p = " << p << " class log likelihood = " << x << " exp -> " << exp(x) << endl;
      double term_a = class_cond / marginal_motif_prior - 1;
      double term_b = exp(-x) * current_class_prior / (1 - marginal_motif_prior);
      double term_c = term_a * term_b;

      // \del \log P(C|X) = P(C) / (P(C|X) * (1 - P(m))) * (P(m|C)/P(m) - 1) * \del P(m|X)

      if(not task.targets.transition.empty())
        t_g[thread_idx] += term_c * t;
      if(not task.targets.emission.empty())
        e_g[thread_idx] += term_c * e;

      if(task.measure == Measure::ClassificationLikelihood) {
        if(not task.targets.transition.empty()) {
          // Compute log likelihood gradients for the full model w.r.t. transition probability
          t_g[thread_idx] += transition_gradient(res.T, task.targets.transition);
        }
        if(not task.targets.emission.empty()) {
          // Compute log likelihood gradients for the full model w.r.t. emission probability
          e_g[thread_idx] += emission_gradient(res.E, task.targets.emission);
        }
       x += res.log_likelihood;
      }
      if(not isfinite(x)) {
        cout << "Error in class likelihood gradient calculation; x is not finite: x = " << x << endl;
        exit(-1);
      }
      l += x;
    }

#pragma omp single
    {
      if(not task.targets.transition.empty())
        for(auto &t: t_g)
          g.transition += t;
      if(not task.targets.emission.empty())
        for(auto &e: e_g)
          g.emission += e;
    }
  }
  if(verbosity >= Verbosity::debug)
    cout << "Data::Set " << dataset.path << " l = " << l << endl;

  return(l);
  */
}

double HMM::site_frequency_difference_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::site_frequenc_difference_gradient(Data::Contrast, Feature)" << endl;
  size_t n_samples = contrast.sets.size();

  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  double signal_counts = 0, control_counts = 0;
  size_t n_signal = 0, n_control = 0;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < n_samples; set_idx++) {
    matrix_t trans_g = zero_matrix(n_states, n_states);
    matrix_t emission_g = zero_matrix(n_states, n_emissions);
    double counts = posterior_gradient(contrast.sets[set_idx], task, present, trans_g, emission_g).posterior;
    bool signal = is_present(contrast.sets[set_idx], present);
    size_t total = contrast.sets[set_idx].set_size + 2 * pseudo_count;
    total += 2 * pseudo_count; // Add pseudo-count
    counts += pseudo_count; // Add pseudo-count
    if(signal) {
      n_signal += total;
      signal_counts += counts;
    } else {
      n_control += total;
      control_counts += counts;
    }

    if(verbosity >= Verbosity::debug)
      cout << "set_idx = " << set_idx << " signal = " << signal << " counts = " << counts << endl;

    if(not task.targets.transition.empty())
      g.transition += (signal ? 1 : -1) * trans_g;
    if(not task.targets.emission.empty())
      g.emission += (signal ? 1 : -1) * emission_g;
  }

  if(verbosity >= Verbosity::verbose) {
    if(not task.targets.transition.empty())
      cout << "Summary posterior transition gradient = " << g.transition << endl;
    if(not task.targets.emission.empty())
      cout << "Summary posterior emission gradient = " << g.emission << endl;
  }

  double delta_freq = signal_counts / n_signal - control_counts / n_control;
  return(delta_freq);
}

double HMM::mutual_information_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information_gradient(Data::Contrast, Feature)" << endl;

  // initialize transition and emission gradient matrices
  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  size_t n_samples = contrast.sets.size();

  // a matrix of expected counts of occurrences and non-occurrences across the samples
  matrix_t counts(n_samples, 2);
  // a vectors of transition and emission gradient matrices for each of the samples
  vector<matrix_t> trans_g(n_samples), emission_g(n_samples);

  // for each of the samples
  for(size_t sample_idx = 0; sample_idx < n_samples; sample_idx++) {
    if(verbosity >= Verbosity::debug)
      cout << "Computing posterior gradient for sample " << contrast.sets[sample_idx].path << endl;

    // compute posterior gradients for transition and emission probabilities, and store the expected occurrences
    counts(sample_idx, 0) = posterior_gradient(contrast.sets[sample_idx], task, present, trans_g[sample_idx], emission_g[sample_idx]).posterior;
    counts(sample_idx, 1) = contrast.sets[sample_idx].set_size - counts(sample_idx, 0);

    if(verbosity >= Verbosity::debug)
      cout << "Transition gradient of sample " << contrast.sets[sample_idx].path << " = " << trans_g[sample_idx] << endl
        << "Emission gradient of sample " << contrast.sets[sample_idx].path << " = " << emission_g[sample_idx] << endl;
  }

  if(verbosity >= Verbosity::debug)
    cout << "Posterior = " << counts << endl;

  // Add pseudo-count
  counts = counts + pseudo_count;

  if(verbosity >= Verbosity::debug)
    cout << "Posterior plus pseudo-count = " << counts << endl;

  double total_r = 0;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++) {
    total_r += counts(set_idx, 1);
    double log_r = log(counts(set_idx, 1));

    if(verbosity >= Verbosity::debug)
      cout << "log_r = " << log_r << endl;

    double factor = log(counts(set_idx, 0)) - log_r;

    if(verbosity >= Verbosity::debug)
      cout << "set_idx = " << set_idx << " factor = " << factor << endl;
    if(not task.targets.transition.empty())
      g.transition += factor * trans_g[set_idx];
    if(not task.targets.emission.empty())
      g.emission += factor * emission_g[set_idx];
  }

  double log_total_r = log(total_r);

  double sum_post = 0;
  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++)
    sum_post += counts(set_idx, 0);

  double factor = log(sum_post) - log_total_r;

  // for each of the samples
  for(size_t set_idx = 0; set_idx < counts.size1(); set_idx++) {
    if(not task.targets.transition.empty())
      g.transition -= factor * trans_g[set_idx];
    if(not task.targets.emission.empty())
      g.emission -= factor * emission_g[set_idx];
  }

  const bool divide_by_data_size = true; // TODO note switch back to 'false' for standard approximate gradient ascent optimization
  if(divide_by_data_size) {
    double n = contrast.set_size + pseudo_count * counts.size1() * counts.size2();
    if(not task.targets.transition.empty())
      g.transition /= n;
    if(not task.targets.emission.empty())
      g.emission /= n;
  }

  if(not task.targets.transition.empty())
    g.transition /= log(2.0);
  if(not task.targets.emission.empty())
    g.emission /= log(2.0);

  if(verbosity >= Verbosity::verbose) {
    if(not task.targets.transition.empty())
      cout << "Summary posterior transition gradient = " << g.transition << endl;
    if(not task.targets.emission.empty())
      cout << "Summary posterior emission gradient = " << g.emission << endl;
  }

  double mi = calc_mutual_information(counts, 0, true, false, false);
//  if(not check_enrichment(contrast, counts, group_idx))
//    mi = -mi;
  return(mi);
}

bool HMM::check_enrichment(const Data::Contrast &contrast, const matrix_t &counts, size_t group_idx) const
{
  string motif = groups[group_idx].name;
  double signal = 0, control = 0;
  double total_signal = 0, total_control = 0;
  for(size_t i = 0; i < counts.size1(); i++) {
    if(contrast.sets[i].motifs.find(motif) != contrast.sets[i].motifs.end()) {
      signal += counts(i,0);
      total_signal += counts(i,0) + counts(i,1);
    } else {
      control += counts(i,0);
      total_control += counts(i,0) + counts(i,1);
    }
  }
  bool ok = signal / total_signal > control / total_control;
  return(ok);
}

double HMM::rank_information_gradient(const Data::Contrast &contrast, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information_gradient(Data::Contrast, Feature)" << endl;

  // initialize transition and emission gradient matrices
  if(not task.targets.transition.empty() and g.transition.size1() == 0)
    g.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and g.emission.size1() == 0)
    g.emission = zero_matrix(n_states, n_emissions);

  double ri = 0;
  for(auto &dataset: contrast)
    ri += rank_information_gradient(dataset, task, present, g);
  return(ri);
}


double HMM::rank_information_gradient(const Data::Set &dataset, const Training::Task &task, bitmask_t present, Gradient &g) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information_gradient(Data::Set, Feature)" << endl;

  const size_t n = dataset.set_size;

  // a vector of expected counts of occurrences across the sequences
  vector_t counts(n);
  // a vectors of transition and emission gradient matrices for each of the sequences
  vector<Gradient> gradients(n);

  // for each of the samples
#pragma omp parallel shared(gradients, counts) if(DO_PARALLEL)
  for(size_t seq_idx = 0; seq_idx < n; seq_idx++) {
    if(verbosity >= Verbosity::debug)
      cout << "Computing posterior gradient for sequence " << dataset.sequences[seq_idx].definition << endl;

    Gradient current_gradient;
    // compute posterior gradients for transition and emission probabilities, and store the expected occurrences
    double current_counts = posterior_gradient(dataset.sequences[seq_idx], task, present, current_gradient.transition, current_gradient.emission).posterior;

    if(verbosity >= Verbosity::debug)
      cout << "Transition gradient of sequence " << dataset.sequences[seq_idx].definition << " = " << gradients[seq_idx].transition << endl
        << "Emission gradient of sequence " << dataset.sequences[seq_idx].definition << " = " << gradients[seq_idx].emission << endl;
#pragma omp critical (store_results)
    {
      if(not task.targets.transition.empty())
        gradients[seq_idx].transition = current_gradient.transition;
      if(not task.targets.emission.empty())
        gradients[seq_idx].emission = current_gradient.emission;
      counts(seq_idx) = current_counts;

    }
  }

  if(verbosity >= Verbosity::debug)
    cout << "Posterior = " << counts << endl;

  Gradient cumul_gradient;
  if(not task.targets.transition.empty())
    cumul_gradient.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty())
    cumul_gradient.emission = zero_matrix(n_states, n_emissions);

  double total = 0;
  for(auto &c: counts)
    total += c;

  if(verbosity >= Verbosity::debug)
    cout << "debug1 " << dataset.path << " " << total << endl;

  double total_factor = 0;
  double cum_count = 0;
  for(size_t i = 0; i < n; i++) {
    if(not task.targets.transition.empty())
      cumul_gradient.transition += gradients[i].transition;
    if(not task.targets.emission.empty())
      cumul_gradient.emission += gradients[i].emission;

    if(i < n-1) {
      cum_count += counts(i);
      double current_factor =
        log(cum_count)
        - log(i + 1 - cum_count)
        + log(n - i - 1 - total + cum_count)
        - log(total - cum_count);
      double current_contribution =
        log(n - total)
        - log(total)
        + log(total - cum_count)
        - log(n - i - 1 - total + cum_count);
      total_factor += current_contribution;

      if(verbosity >= Verbosity::debug)
        cout << "debug2 " << dataset.path << " " << i << " " << cum_count << " " << current_factor << " " << current_contribution << " " << total_factor << endl;

      if(not task.targets.transition.empty())
        g.transition += cumul_gradient.transition * current_factor;
      if(not task.targets.emission.empty())
        g.emission += cumul_gradient.emission * current_factor;
    }
  }


  if(not task.targets.transition.empty())
    g.transition += cumul_gradient.transition * total_factor;
  if(not task.targets.emission.empty())
    g.emission += cumul_gradient.emission * total_factor;

  if(not task.targets.transition.empty())
    g.transition /=  log(2.0) * n * n; // we compute the mean mutual information over all ranks, rather than the sum
    // g.transition /=  log(2.0) * n;
  if(not task.targets.emission.empty())
    g.emission /=  log(2.0) * n * n; // we compute the mean mutual information over all ranks, rather than the sum
    // g.emission /=  log(2.0) * n;

    // TODO: reactivate
    // Add pseudo-count
    // counts = counts + pseudo_count;

    // if(verbosity >= Verbosity::debug)
    //   cout << "Posterior plus pseudo-count = " << counts << endl;
  if(verbosity >= Verbosity::debug) {
    if(not task.targets.transition.empty())
      cout << "Summary rank mutual information transition gradient = " << g.transition << endl;
    if(not task.targets.emission.empty())
      cout << "Summary rank mutual information emission gradient = " << g.emission << endl;
  }

  //  cout << "TODO: implement." << endl;
  //  exit(-1);
  double ri = calc_rank_information(counts, pseudo_count);
  return(ri);
}

/** This computes the log likelihood gradient w.r.t. the transformed transition probabilities.*/
matrix_t HMM::transition_gradient(const matrix_t &T, const Training::Range &range) const
{
  matrix_t m = zero_matrix(n_states, n_states);
  for(auto i: range)
    for(auto j: succ[i])
      for(auto k: succ[i])
        m(i,j) += T(i,k) * (((j==k) ? 1 : 0) - transition(i,j));
  return(m);
}

/** This computes the log likelihood gradient w.r.t. the transformed emission probabilities.*/
matrix_t HMM::emission_gradient(const matrix_t &E, const Training::Range &range) const
{
  matrix_t m = zero_matrix(n_states, n_emissions);
  for(auto j: range)
    for(size_t k = 0; k < n_emissions; k++)
      for(size_t l = 0; l < n_emissions; l++)
        m(j,k) += E(j,l) * (((k==l) ? 1 : 0) - emission(j,k));
  return(m);
}

HMM::posterior_gradient_t HMM::posterior_gradient(const Data::Seq &seq, const Training::Task &task, bitmask_t present, matrix_t &transition_g, matrix_t &emission_g) const
{
  // Training::Task task(task_);
  if(verbosity >= Verbosity::debug)
    cout << "Posterior gradient calculation (Seq, Feature)." << endl;

  SubHMM subhmm(*this, complementary_states_mask(present));

  if(verbosity >= Verbosity::debug)
    cout << *this << endl
      << subhmm << endl;

  if(verbosity >= Verbosity::debug) {
    cout << "Transition targets are";
    for(auto x:task.targets.transition)
      cout << " " << x;
    cout << endl;
    cout << "Emission targets are";
    for(auto x:task.targets.emission)
      cout << " " << x;
    cout << endl;
  }

  if(not task.targets.transition.empty())
    transition_g = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty())
    emission_g = zero_matrix(n_states, n_emissions);

  Training::Targets reduced_targets = subhmm.map_down(task.targets);

  if(verbosity >= Verbosity::debug) {
    cout << "targets emission = ";
    for(auto &x: task.targets.emission)
      cout << " " << x;
    cout << endl;
    cout << "targets transition = ";
    for(auto &x: task.targets.transition)
      cout << " " << x;
    cout << endl;
    cout << "reduced targets emission = ";
    for(auto &x: reduced_targets.emission)
      cout << " " << x;
    cout << endl;
    cout << "reduced targets transition = ";
    for(auto &x: reduced_targets.transition)
      cout << " " << x;
    cout << endl;
  }

  // Compute gradient for each sequence
//  int thread_idx = omp_get_thread_num();
//  if(verbosity >= Verbosity::debug)
//    cout << "Thread " << thread_idx << " Data sample " <<  endl << seq2string(seq.sequence) << endl;

  // Compute expected statistics, for the full and reduced models
  matrix_t T, Tr, E, Er;
  double logp = BaumWelchIteration_single(T, E, seq, task.targets);
  double logpr = subhmm.BaumWelchIteration_single(Tr, Er, seq, reduced_targets);

  Tr = subhmm.lift_transition(Tr);
  Er = subhmm.lift_emission(Er);

  if(verbosity >= Verbosity::debug)
    cout << "Full logp = " << logp << endl
      << "Reduced logp = " << logpr << endl
      << "Expected transitions full = " << T << endl
      << "Expected emissions full = " << E << endl
      << "Expected transitions constitutive_range = " << Tr << endl
      << "Expected emissions constitutive_range = " << Er << endl;

  if(not task.targets.transition.empty()) {
    // Compute log likelihood gradients for the full model w.r.t. transition probability 
    matrix_t t = transition_gradient(T, task.targets.transition);
    // Compute log likelihood gradients for the reduced model w.r.t. transition probability
    matrix_t tr = transition_gradient(Tr, task.targets.transition);

    // Compute posterior probability gradients for the reduced model w.r.t. transition probability and accumulate
    transition_g += exp(logpr - logp) * (t - tr);
  }

  if(not task.targets.emission.empty()) {
    // Compute log likelihood gradients for the full model w.r.t. emission probability
    matrix_t e = emission_gradient(E, task.targets.emission);
    // Compute log likelihood gradients for the reduced model w.r.t. emission probability
    matrix_t er = emission_gradient(Er, task.targets.emission);

    // Compute posterior probability gradients for the reduced model w.r.t. emission probability and accumulate
    emission_g += exp(logpr - logp) * (e - er);
  }

  double posterior = 1 - exp(logpr - logp);

  if(verbosity >= Verbosity::verbose)
    cout << "The posterior coming from the gradient calculus: " << posterior << endl;
  posterior_gradient_t result = {logp, posterior, T, E};
  return(result);
}

HMM::posterior_t HMM::posterior_gradient(const Data::Set &dataset, const Training::Task &task, bitmask_t present, matrix_t &transition_g, matrix_t &emission_g) const
{
  // Training::Task task(task_);
  if(verbosity >= Verbosity::verbose)
    cout << "Posterior gradient calculation (Feature)." << endl;

  SubHMM subhmm(*this, complementary_states_mask(present));

  // cout << subhmm << endl;

  if(verbosity >= Verbosity::debug) {
    cout << "Transition targets are";
    for(auto x:task.targets.transition)
      cout << " " << x;
    cout << endl;
    cout << "Emission targets are";
    for(auto x:task.targets.emission)
      cout << " " << x;
    cout << endl;
  }

  if(not task.targets.transition.empty())
    transition_g = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty())
    emission_g = zero_matrix(n_states, n_emissions);
  double posterior = 0;
  vector<matrix_t> t_g, e_g;

  Training::Targets reduced_targets = subhmm.map_down(task.targets);
 
  if(verbosity >= Verbosity::debug) {
    cout << "targets emission = ";
    for(auto &x: task.targets.emission)
      cout << " " << x;
    cout << endl;
    cout << "targets transition = ";
    for(auto &x: task.targets.transition)
      cout << " " << x;
    cout << endl;
    cout << "reduced targets emission = ";
    for(auto &x: reduced_targets.emission)
      cout << " " << x;
    cout << endl;
    cout << "reduced targets transition = ";
    for(auto &x: reduced_targets.transition)
      cout << " " << x;
    cout << endl;
  }

  double l = 0;

#pragma omp parallel shared(emission_g, transition_g) if(DO_PARALLEL)
  {
    // Initalize storage for thread intermediate results
#pragma omp single
    {
      size_t n_threads = omp_get_num_threads();
      //    cout << "B Num threads = " << n_threads << endl;
      if(not task.targets.transition.empty())
        t_g = vector<matrix_t> (n_threads, zero_matrix(transition_g.size1(), transition_g.size2()));
      if(not task.targets.emission.empty())
        e_g = vector<matrix_t> (n_threads, zero_matrix(emission_g.size1(), emission_g.size2()));
    }

    // Compute gradient for each sequence
#pragma omp for reduction (+:posterior, l)
    for(size_t i = 0; i < dataset.sequences.size(); i++) {
      int thread_idx = omp_get_thread_num();
      if(verbosity >= Verbosity::debug)
        cout << "Thread " << thread_idx << " Data sample " << i << endl << seq2string(dataset.sequences[i].isequence) << endl;

      // Compute expected statistics, for the full and reduced models
      matrix_t T, Tr, E, Er;
      double logp = BaumWelchIteration_single(T, E, dataset.sequences[i], task.targets);
      double logpr = subhmm.BaumWelchIteration_single(Tr, Er, dataset.sequences[i], reduced_targets);

      Tr = subhmm.lift_transition(Tr);
      Er = subhmm.lift_emission(Er);

      if(verbosity >= Verbosity::debug)
        cout << "Full logp = " << logp << endl
          << "Reduced logp = " << logpr << endl
          << "Expected transitions full = " << T << endl
          << "Expected emissions full = " << E << endl
          << "Expected transitions constitutive_range = " << Tr << endl
          << "Expected emissions constitutive_range = " << Er << endl;

      if(not task.targets.transition.empty()) {
        // Compute log likelihood gradients for the full model w.r.t. transition probability 
        matrix_t t = transition_gradient(T, task.targets.transition);
        // Compute log likelihood gradients for the reduced model w.r.t. transition probability
        matrix_t tr = transition_gradient(Tr, task.targets.transition);

        // Compute posterior probability gradients for the reduced model w.r.t. transition probability and accumulate
        t_g[thread_idx] += exp(logpr - logp) * (t - tr);
      }

      if(not task.targets.emission.empty()) {
        // Compute log likelihood gradients for the full model w.r.t. emission probability
        matrix_t e = emission_gradient(E, task.targets.emission);
        // Compute log likelihood gradients for the reduced model w.r.t. emission probability
        matrix_t er = emission_gradient(Er, task.targets.emission);

        // Compute posterior probability gradients for the reduced model w.r.t. emission probability and accumulate
        e_g[thread_idx] += exp(logpr - logp) * (e - er);
      }

      posterior += 1 - exp(logpr - logp);
      l += logp;
    }

    // Collect results of threads
#pragma omp single
    {
      if(not task.targets.transition.empty())
        for(auto &x: t_g)
          transition_g += x;
      if(not task.targets.emission.empty())
        for(auto &x: e_g)
          emission_g += x;
    }
  }

  if(verbosity >= Verbosity::verbose)
    cout << "The posterior coming from the gradient calculus: " << posterior << endl;
  posterior_t result = {l, posterior};
  return(result);
}

