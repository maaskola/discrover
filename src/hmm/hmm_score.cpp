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
 *       Filename:  hmm_score.cpp
 *
 *    Description:  HMM routines that calculate various statistical scores.
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include "../aux.hpp"
#include "hmm.hpp"
#include "subhmm.hpp"

using namespace std;

#define DO_PARALLEL 1

size_t make_mask(const vector<size_t> &v) {
  size_t x = 0;
  for(auto y: v)
    x |= 1 << y;
  return(x);
}

double HMM::mutual_information(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing feature-wise mutual information." << endl;
  size_t present_mask = make_mask(present_groups); size_t absent_mask = make_mask(absent_groups);
  return(mutual_information(data, present_mask, absent_mask));
}

double HMM::rank_information(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing feature-wise mutual information." << endl;
  size_t present_mask = make_mask(present_groups); size_t absent_mask = make_mask(absent_groups);
  return(rank_information(data, present_mask, absent_mask));
}

double HMM::matthews_correlation_coefficient(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing sum of feature-wise Matthew's correlation coefficient." << endl;
  size_t present_mask = make_mask(present_groups); size_t absent_mask = make_mask(absent_groups);
  return(matthews_correlation_coefficient(data, present_mask, absent_mask));
}

// TODO FIX ABSENT - done?
confusion_matrix reduce(const vector_t &v, size_t present_mask, const Data::Series &data, const vector<Group> &groups, bool word_stats) {
  confusion_matrix m = {0, 0, 0, 0};
  for(size_t sample_idx = 0; sample_idx < v.size(); sample_idx++) {
    bool signal = false;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if(((1 << group_idx) & present_mask) != 0 and
          data.sets[sample_idx].motifs.find(groups[group_idx].name) != end(data.sets[sample_idx].motifs)) {
        signal = true;
        break;
      }
    if(signal) {
      m.true_positives += v[sample_idx];
      m.false_negatives += (word_stats ? data.sets[sample_idx].seq_size : data.sets[sample_idx].set_size) - v[sample_idx];
    } else {
      m.false_positives += v[sample_idx];
      m.true_negatives += (word_stats ? data.sets[sample_idx].seq_size : data.sets[sample_idx].set_size) - v[sample_idx];
    }
  }
  return(m);
}

double HMM::mutual_information(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information(Data::Series, Feature)" << endl;
  vector_t posterior = posterior_atleast_one(data, present_mask, absent_mask);
  matrix_t m(posterior.size(), 2);
  for(size_t i = 0; i < posterior.size(); i++) {
    m(i,0) = posterior(i);
    m(i,1) = data.sets[i].set_size - posterior(i);
  }
  m = m + pseudo_count;
  cout << "HMM::mutual_information(Data::Series, Feature) present = " << present_mask << " absent = " << absent_mask << endl
    << "counts = " << m << endl;
  double mi = calc_mutual_information(m, 0, true, false, false);
//  if(not check_enrichment(data, m, group_idx))
//    mi = -mi;
  return(mi);
}

double HMM::rank_information(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Series, Feature)" << endl;
  double ri = 0;
  for(auto &set: data)
    ri += rank_information(set, present_mask, absent_mask);
  return(ri);
}


double HMM::rank_information(const Data::Set &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Set, Feature)" << endl;
  vector_t posterior = posterior_atleast_one(data, present_mask, absent_mask);
  return(calc_rank_information(posterior, pseudo_count));
}

double HMM::matthews_correlation_coefficient(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  // TODO find out if there's a proper generalization of the MCC to multiple experiments
  vector_t posterior = posterior_atleast_one(data, present_mask, absent_mask);
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present_mask, data, groups, false) + pseudo_count;
  return(calc_matthews_correlation_coefficient(m));
}

double HMM::log_likelihood(const Data::Collection &data) const
{
  double l = 0;
  for(auto &series: data)
    l += log_likelihood(series);
  return(l);
}

double HMM::log_likelihood(const Data::Series &data) const
{
  double l = 0;
  for(auto &data_set: data)
    l += log_likelihood(data_set);
  return(l);
}

double HMM::log_likelihood(const Data::Set &data) const
{
  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++)
    l += log_likelihood_from_scale(compute_forward_scale(data.sequences[i]));
  return(l);
}

vector_t HMM::expected_posterior(const Data::Series &data, size_t group_idx) const
{
  vector_t v = zero_vector(data.sets.size());
  for(size_t i = 0; i < data.sets.size(); i++)
    v(i) = expected_posterior(data.sets[i], group_idx);
  return(v);
};

double HMM::expected_posterior(const Data::Set &data, size_t group_idx) const
{
  double m = 0;
#pragma omp parallel for schedule(static) reduction(+:m) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++) {
    vector_t scale;
    matrix_t f = compute_forward_scaled(data.sequences[i], scale);
    matrix_t b = compute_backward_prescaled(data.sequences[i], scale);
    m += expected_state_posterior(groups[group_idx].states[0], f, b, scale); // Assume the first state of each motif is constitutive for the motif
  }
  return(m);
};

double HMM::expected_posterior(const Data::Seq &data, size_t group_idx) const
{
  vector_t scale;
  matrix_t f = compute_forward_scaled(data, scale);
  matrix_t b = compute_backward_prescaled(data, scale);
  double m = expected_state_posterior(groups[group_idx].states[0], f, b, scale); // Assume the first state of each motif is constitutive for the motif
  return(m);
};

vector_t HMM::posterior_atleast_one(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug) {
    cout << "HMM::posterior_atleast_one(Data::Series, present_groups=(";
    bool first = true; for(auto &x: present_groups) { cout << (first ? "" : ",") << x; first = false; } cout << ") absent_groups=(";
    first = true; for(auto &x: absent_groups) { cout << (first ? "" : ",") << x; first = false; } cout << ")" << endl;
  }
  return(posterior_atleast_one(data, make_mask(present_groups), make_mask(absent_groups)));
}

vector_t HMM::posterior_atleast_one(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Series, present_mask=" << present_mask << ", absent_mask=" << absent_mask << endl;
  vector_t v(data.sets.size());
  for(size_t i = 0; i < data.sets.size(); i++)
    v[i] = sum_posterior_atleast_one(data.sets[i], present_mask, absent_mask);
  return(v);
}

vector_t HMM::posterior_atleast_one(const Data::Set &data, size_t present_mask, size_t absent_mask) const
{
  // if(verbosity >= Verbosity::debug)
  if(true)
    cout << "HMM::posterior_atleast_one(Data::Set = " << data.path << ", present_mask = " << present_mask << ", absent_mask = " << absent_mask << endl;

  cout << "complementary_states_mask(present_mask)) =";
  for(auto &state: complementary_states_mask(present_mask))
    cout << " " << state;
  cout << endl;
  cout << "complementary_states_mask(present_mask | absent_mask)) =";
  for(auto &state: complementary_states_mask(present_mask | absent_mask))
    cout << " " << state;
  cout << endl;

  // TODO FIX ABSENT - done?
  vector_t vec(data.set_size);
  if(absent_mask == 0) {
    SubHMM subhmm(*this, complementary_states_mask(present_mask));
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
    for(size_t i = 0; i < data.set_size; i++) {
      double logp = log_likelihood_from_scale(compute_forward_scale(data.sequences[i]));
      double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(data.sequences[i]));
      double z = 1 - exp(logp_wo_motif - logp);
      if(verbosity >= Verbosity::debug)
        cout << "seq = " << data.sequences[i].definition
          /*  << " " << data.sequences[i].sequence */
          << " logp = " << logp
          << " logp_wo_motif = " << logp_wo_motif
          << " z = " << z << endl;
      vec[i] = z;
    }
  } else {
    SubHMM subhmm_wo_abs(*this, complementary_states_mask(absent_mask));
    SubHMM subhmm_wo_abs_wo_motif(*this, complementary_states_mask(present_mask | absent_mask));
    cout << "Full     :" << *this << endl;
    cout << "Reduced 1:" << subhmm_wo_abs << endl;
    cout << "Reduced 2:" << subhmm_wo_abs_wo_motif << endl;
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
    for(size_t i = 0; i < data.set_size; i++) {
      double logp = log_likelihood_from_scale(compute_forward_scale(data.sequences[i]));
      double logp_wo_abs = log_likelihood_from_scale(subhmm_wo_abs.compute_forward_scale(data.sequences[i]));
      double logp_wo_abs_wo_motif = log_likelihood_from_scale(subhmm_wo_abs_wo_motif.compute_forward_scale(data.sequences[i]));

      double z = exp(logp_wo_abs - logp) - exp(logp_wo_abs_wo_motif - logp);
      if(verbosity >= Verbosity::debug) {
      // if(true) {
        stringstream s;
        s << "seq = " << data.sequences[i].definition
          /*  << " " << data.sequences[i].sequence */
          << " logp = " << logp
          << " logp_wo_abs = " << logp_wo_abs
          << " logp_wo_abs_wo_motif = " << logp_wo_abs_wo_motif
          << " z = " << z << endl;
        cout << s.str();
      }
      vec[i] = z;
    }
  }

  // if(verbosity >= Verbosity::debug)
  if(true)
    cout << "HMM::posterior_atleast_one(Data::Set = " << data.path << ", present_mask = " << present_mask << ", absent_mask = " << absent_mask << " vec = " << vec << endl;
  return(vec);
};



double HMM::sum_posterior_atleast_one(const Data::Set &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_posterior_atleast_one(Data::Set = " << data.path << ", present_mask = " << present_mask << ", absent_mask = " << absent_mask << endl;

  vector_t counts = posterior_atleast_one(data, present_mask, absent_mask);

  double m = 0;
  for(auto &x: counts)
    m += x;

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Set = " << data.path << ", present_mask = " << present_mask << ", absent_mask = " << absent_mask << " m = " << m << endl;
  return(m);
};

HMM::posterior_t HMM::posterior_atleast_one(const Data::Seq &data, size_t present_mask, size_t absent_mask) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq, present_mask = " << present_mask << ", absent_mask = " << absent_mask << endl;

  // TODO FIX ABSENT - done?

  double logp = log_likelihood_from_scale(compute_forward_scale(data));

  double z;
  if(absent_mask == 0) {
    SubHMM subhmm(*this, complementary_states_mask(present_mask));
    double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(data));

    z = 1 - exp(logp_wo_motif - logp);
    if(verbosity >= Verbosity::debug)
      cout << "seq = " << data.definition << " " << data.sequence << " logp = " << logp << " logp_wo_motif = " << logp_wo_motif << " z = " << z << endl;
  } else {
    SubHMM subhmm_wo_abs(*this, complementary_states_mask(absent_mask));
    SubHMM subhmm_wo_abs_wo_motif(*this, complementary_states_mask(present_mask | absent_mask));
    if(verbosity >= Verbosity::debug) {
      cout << "Full     :" << *this << endl;
      cout << "Reduced 1:" << subhmm_wo_abs << endl;
      cout << "Reduced 2:" << subhmm_wo_abs_wo_motif << endl;
    }
    double logp_wo_abs = log_likelihood_from_scale(subhmm_wo_abs.compute_forward_scale(data));
    double logp_wo_abs_wo_motif = log_likelihood_from_scale(subhmm_wo_abs_wo_motif.compute_forward_scale(data));

    z = exp(logp_wo_abs - logp) - exp(logp_wo_abs_wo_motif - logp);
    if(verbosity >= Verbosity::debug) {
      stringstream s;
      s << "seq = " << data.definition
        /*  << " " << data.sequences[i].sequence */
        << " logp = " << logp
        << " logp_wo_abs = " << logp_wo_abs
        << " logp_wo_abs_wo_motif = " << logp_wo_abs_wo_motif
        << " z = " << z << endl;
      cout << s.str();
    }
  }

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq, present_mask = " << present_mask << ", absent_mask = " << absent_mask << " z = " << z << endl;
  posterior_t res = {logp, z};
  return(res);
};

vector_t HMM::viterbi_atleast_one(const Data::Series &data, size_t group_idx) const
{
  vector_t v(data.sets.size());
  for(size_t i = 0; i < data.sets.size(); i++)
    v[i] = viterbi_atleast_one(data.sets[i], group_idx);
  return(v);
}

double HMM::viterbi_atleast_one(const Data::Set &data, size_t group_idx) const
{
  double m = 0;
#pragma omp parallel for schedule(static) reduction(+:m) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++) {
    StatePath path;
    viterbi(data.sequences[i], path);
    if(find(path.begin(), path.end(), groups[group_idx].states[0]) != path.end())
      m++;
  }
  return(m);
};

double HMM::compute_score_all_motifs(const Data::Collection &data, const Measures::Continuous::Measure &measure, bool weighting) const
{
  vector<size_t> all_motifs;
  for(size_t i = 0; i < groups.size(); i++)
    if(groups[i].kind == Group::Kind::Motif)
      all_motifs.push_back(i);
  return(compute_score(data, measure, weighting, all_motifs, vector<size_t>()));
}

double HMM::compute_score(const Data::Collection &data, const Measures::Continuous::Measure &measure, bool weighting, const vector<size_t> &present_motifs, const vector<size_t> &absent_motifs) const
{
  double score = 0;
  double W = 0;
  double w;
  switch(measure) {
    case Measure::Likelihood:
      score = log_likelihood(data);
      break;
    case Measure::MutualInformation:
      for(auto &series: data) {
        W += w = series.set_size;
        score += mutual_information(series, present_motifs, absent_motifs) * (weighting ? w : 1);
      }
      if(weighting)
        score /= W;
      break;
    case Measure::RankInformation:
      for(auto &series: data)
        score = rank_information(series, present_motifs, absent_motifs);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      for(auto &series: data)
        score += matthews_correlation_coefficient(series, present_motifs, absent_motifs);
      break;
    case Measure::DeltaFrequency:
      for(auto &series: data)
        score += dips_sitescore(series, present_motifs, absent_motifs);
      break;
    case Measure::LogLikelihoodDifference:
      for(auto &series: data)
        score += log_likelihood_difference(series, present_motifs, absent_motifs);
      break;
    case Measure::ClassificationPosterior:
      for(auto &series: data)
        score += class_likelihood(series, present_motifs, absent_motifs, true);
      break;
    case Measure::ClassificationLikelihood:
      for(auto &series: data)
        score += class_likelihood(series, present_motifs, absent_motifs, false);
      break;
    default:
      cout << "Score calculation for '" << measure2string(measure) << "' is not implemented." << endl;
      assert(0);
  }
  return(score);
}

double HMM::class_likelihood(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups, bool compute_posterior) const
{
  return(class_likelihood(data, make_mask(present_groups), make_mask(absent_groups), compute_posterior));
}

double HMM::class_likelihood(const Data::Series &data, size_t present_mask, size_t absent_mask, bool compute_posterior) const
{
  double l = 0;
  for(auto &set: data)
    l += class_likelihood(set, present_mask, absent_mask, compute_posterior);
  if(verbosity >= Verbosity::debug)
    cout << "Data::Series l = " << l << endl;
  return(l);
}


double HMM::class_likelihood(const Data::Set &data, size_t present_mask, size_t absent_mask, bool compute_posterior) const
{
  // TODO FIX ABSENT
  // const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
  const double marginal_motif_prior  = 0.5;
  // TODO FIX ABSENT
  // const double class_cond = get_class_motif_prior(data.sha1, group_idx);
  const double class_cond = 0.2;
  const double log_class_prior = log(get_class_prior(data.sha1));

  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++) {
    posterior_t res = posterior_atleast_one(data.sequences[i], present_mask, absent_mask);
    double p = res.posterior;
    double x = log_class_prior + log(p * class_cond / marginal_motif_prior + (1-p) * (1-class_cond) / (1-marginal_motif_prior));
    if(not compute_posterior)
      x += res.log_likelihood;
    if(verbosity >= Verbosity::debug)
      cout << "Sequence " << data.sequences[i].definition << " p = " << p << " class log likelihood = " << x << " exp -> " << exp(x) << endl;
    l += x;
  }
  if(verbosity >= Verbosity::debug)
    cout << "Data::Set " << data.path << " l = " << l << endl;
  return(l);
}

double HMM::log_likelihood_difference(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  return(log_likelihood_difference(data, make_mask(present_groups), make_mask(absent_groups)));
}

double HMM::log_likelihood_difference(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  double d = 0;
  for(size_t sample_idx = 0; sample_idx < data.sets.size(); sample_idx++) {
    // TODO FIX ABSENT - done?
    bool signal = false;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if(((1 << group_idx) & present_mask) != 0 and
          data.sets[sample_idx].motifs.find(groups[group_idx].name) != end(data.sets[sample_idx].motifs)) {
        signal = true;
        break;
      }
    d += (signal ? 1 : -1) * log_likelihood(data.sets[sample_idx]);
  }
  return(d);
}

double HMM::dips_sitescore(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  return(dips_sitescore(data, make_mask(present_groups), make_mask(absent_groups)));
}

double HMM::dips_sitescore(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  vector_t posterior = posterior_atleast_one(data, present_mask, absent_mask);
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present_mask, data, groups, false) + pseudo_count;
  size_t signal_size = m.true_positives + m.false_negatives;
  size_t control_size = m.false_positives + m.true_negatives;

  double t = m.true_positives / signal_size - m.false_positives / control_size;
  return(t);
}

double HMM::dips_tscore(const Data::Series &data, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  return(dips_tscore(data, make_mask(present_groups), make_mask(absent_groups)));
}

double HMM::dips_tscore(const Data::Series &data, size_t present_mask, size_t absent_mask) const
{
  // TODO FIX ABSENT
  // vector_t posterior = expected_posterior(data, group_idx);
  vector_t posterior;
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present_mask, data, groups, true) + pseudo_count;
  size_t signal_size = m.true_positives + m.false_negatives;
  size_t control_size = m.false_positives + m.true_negatives;

  double t = m.true_positives / signal_size - m.false_positives / control_size;
  return(t);
}

double HMM::information_content(size_t motif_idx) const
{
  double h = 0;
  for(auto i: groups[motif_idx].states)
    for(size_t j = 0; j < n_emissions; j++)
      if(emission(i,j) != 0)
        h -= emission(i,j) * log(emission(i,j));
  h /= log(2.0);
  double ic = 2 * groups[motif_idx].states.size() - h;
  return(ic);
}

