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

HMM::bitmask_t make_mask(const vector<size_t> &v) {
  HMM::bitmask_t x = 0;
  for(auto y: v) {
    if(y > HMM::max_motifs) {
      cout << "Error: trying to construct mask for too many motifs! The offending index was: " << y << ", and there may only be " << HMM::max_motifs << " motifs in this version." << endl;
      exit(-1);
    }
    // TODO implement in terms of set() or operator[] methods
    x |= 1 << y;
  }
  return(x);
}

vector<size_t> unpack_mask(const HMM::bitmask_t x) {
  vector<size_t> v;
  size_t z = 0;
  HMM::bitmask_t y = 1;
  while(x.to_ullong() >= y.to_ullong()) {
    if((x & y) != 0)
      v.push_back(z);
    // TODO implement in terms of test() or operator[] methods
    y = y << 1;
    z++;
  }
  return(v);
}

double HMM::mutual_information(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing mutual information." << endl;
  bitmask_t present = make_mask(present_groups); bitmask_t absent = make_mask(absent_groups);
  return(mutual_information(contrast, present, absent));
}

double HMM::rank_information(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing mutual information." << endl;
  bitmask_t present = make_mask(present_groups); bitmask_t absent = make_mask(absent_groups);
  return(rank_information(contrast, present, absent));
}

double HMM::matthews_correlation_coefficient(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing sum of Matthew's correlation coefficient." << endl;
  bitmask_t present = make_mask(present_groups); bitmask_t absent = make_mask(absent_groups);
  return(matthews_correlation_coefficient(contrast, present, absent));
}

// TODO FIX ABSENT - done?
confusion_matrix reduce(const vector_t &v, HMM::bitmask_t present, const Data::Contrast &contrast, const vector<Group> &groups, bool word_stats) {
  confusion_matrix m = {0, 0, 0, 0};
  for(size_t sample_idx = 0; sample_idx < v.size(); sample_idx++) {
    bool signal = false;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if((HMM::bitmask_t(1 << group_idx) & present) != 0 and
          contrast.sets[sample_idx].motifs.find(groups[group_idx].name) != end(contrast.sets[sample_idx].motifs)) {
        signal = true;
        break;
      }
    if(signal) {
      m.true_positives += v[sample_idx];
      m.false_negatives += (word_stats ? contrast.sets[sample_idx].seq_size : contrast.sets[sample_idx].set_size) - v[sample_idx];
    } else {
      m.false_positives += v[sample_idx];
      m.true_negatives += (word_stats ? contrast.sets[sample_idx].seq_size : contrast.sets[sample_idx].set_size) - v[sample_idx];
    }
  }
  return(m);
}

double HMM::mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information(Data::Contrast)" << endl;
  vector_t posterior = posterior_atleast_one(contrast, present, absent);
  matrix_t m(posterior.size(), 2);
  for(size_t i = 0; i < posterior.size(); i++) {
    m(i,0) = posterior(i);
    m(i,1) = contrast.sets[i].set_size - posterior(i);
  }
  m = m + pseudo_count;
  cout << "HMM::mutual_information(Data::Contrast)" << endl
    << "present = " << present << endl
    << "absent  = " << absent << endl
    << "counts  = " << m << endl;
  double mi = calc_mutual_information(m, 0, true, false, false);
//  if(not check_enrichment(contrast, m, group_idx))
//    mi = -mi;
  return(mi);
}

double HMM::rank_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Contrast)" << endl;
  double ri = 0;
  for(auto &dataset: contrast)
    ri += rank_information(dataset, present, absent);
  return(ri);
}


double HMM::rank_information(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Set)" << endl;
  vector_t posterior = posterior_atleast_one(dataset, present, absent);
  return(calc_rank_information(posterior, pseudo_count));
}

double HMM::matthews_correlation_coefficient(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  // TODO find out if there's a proper generalization of the MCC to multiple experiments
  vector_t posterior = posterior_atleast_one(contrast, present, absent);
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present, contrast, groups, false) + pseudo_count;
  return(calc_matthews_correlation_coefficient(m));
}

double HMM::log_likelihood(const Data::Collection &collection) const
{
  double l = 0;
  for(auto &contrast: collection)
    l += log_likelihood(contrast);
  return(l);
}

double HMM::log_likelihood(const Data::Contrast &contrast) const
{
  double l = 0;
  for(auto &dataset: contrast)
    l += log_likelihood(dataset);
  return(l);
}

double HMM::log_likelihood(const Data::Set &dataset) const
{
  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < dataset.set_size; i++)
    l += log_likelihood_from_scale(compute_forward_scale(dataset.sequences[i]));
  return(l);
}

vector_t HMM::expected_posterior(const Data::Contrast &contrast, bitmask_t present) const
{
  vector_t v = zero_vector(contrast.sets.size());
  for(size_t i = 0; i < contrast.sets.size(); i++)
    v(i) = expected_posterior(contrast.sets[i], present);
  return(v);
};

double HMM::expected_posterior(const Data::Set &dataset, bitmask_t present) const
{
  double m = 0;
  vector<size_t> present_groups = unpack_mask(present);
#pragma omp parallel for schedule(static) reduction(+:m) if(DO_PARALLEL)
  for(size_t i = 0; i < dataset.set_size; i++) {
    vector_t scale;
    matrix_t f = compute_forward_scaled(dataset.sequences[i], scale);
    matrix_t b = compute_backward_prescaled(dataset.sequences[i], scale);
    for(auto group_idx: present_groups)
      m += expected_state_posterior(groups[group_idx].states[0], f, b, scale); // Assume the first state of each motif is constitutive for the motif
  }
  return(m);
};

double HMM::expected_posterior(const Data::Seq &seq, bitmask_t present) const
{
  vector<size_t> present_groups = unpack_mask(present);
  vector_t scale;
  matrix_t f = compute_forward_scaled(seq, scale);
  matrix_t b = compute_backward_prescaled(seq, scale);
  double m = 0;
  for(auto group_idx: present_groups)
    m += expected_state_posterior(groups[group_idx].states[0], f, b, scale); // Assume the first state of each motif is constitutive for the motif
  return(m);
};

vector_t HMM::posterior_atleast_one(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Contrast)" << endl
      << "present =" << present << endl
      << "absent  =" << absent << endl;
  vector_t v(contrast.sets.size());
  for(size_t i = 0; i < contrast.sets.size(); i++)
    v[i] = sum_posterior_atleast_one(contrast.sets[i], present, absent);
  return(v);
}

vector_t HMM::posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "absent  = " << absent << endl;

  if(verbosity >= Verbosity::debug) {
    cout << "complementary_states_mask(present)) =";
    for(auto &state: complementary_states_mask(present))
      cout << " " << state;
    cout << endl;
    cout << "complementary_states_mask(present | absent)) =";
    for(auto &state: complementary_states_mask(present | absent))
      cout << " " << state;
    cout << endl;
  }

  // TODO FIX ABSENT - done?
  vector_t vec(dataset.set_size);
  if(absent == 0) {
    SubHMM subhmm(*this, complementary_states_mask(present));
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
    for(size_t i = 0; i < dataset.set_size; i++) {
      double logp = log_likelihood_from_scale(compute_forward_scale(dataset.sequences[i]));
      double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(dataset.sequences[i]));
      double z = 1 - exp(logp_wo_motif - logp);
      if(verbosity >= Verbosity::debug)
        cout << "seq = " << dataset.sequences[i].definition
          /*  << " " << dataset.sequences[i].sequence */
          << " logp = " << logp
          << " logp_wo_motif = " << logp_wo_motif
          << " z = " << z << endl;
      vec[i] = z;
    }
  } else {
    SubHMM subhmm_wo_abs(*this, complementary_states_mask(absent));
    SubHMM subhmm_wo_abs_wo_motif(*this, complementary_states_mask(present | absent));
    if(verbosity >= Verbosity::debug) {
      cout << "Full     :" << *this << endl;
      cout << "Reduced 1:" << subhmm_wo_abs << endl;
      cout << "Reduced 2:" << subhmm_wo_abs_wo_motif << endl;
    }
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
    for(size_t i = 0; i < dataset.set_size; i++) {
      double logp = log_likelihood_from_scale(compute_forward_scale(dataset.sequences[i]));
      double logp_wo_abs = log_likelihood_from_scale(subhmm_wo_abs.compute_forward_scale(dataset.sequences[i]));
      double logp_wo_abs_wo_motif = log_likelihood_from_scale(subhmm_wo_abs_wo_motif.compute_forward_scale(dataset.sequences[i]));

      double z = exp(logp_wo_abs - logp) - exp(logp_wo_abs_wo_motif - logp);
      if(verbosity >= Verbosity::debug) {
        stringstream s;
        s << "seq = " << dataset.sequences[i].definition
          /*  << " " << dataset.sequences[i].sequence */
          << " logp = " << logp
          << " logp_wo_abs = " << logp_wo_abs
          << " logp_wo_abs_wo_motif = " << logp_wo_abs_wo_motif
          << " z = " << z << endl;
        cout << s.str();
      }
      vec[i] = z;
    }
  }

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "absent  = " << absent << endl
      << "vec = " << vec << endl;
  return(vec);
};



double HMM::sum_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "absent  = " << absent << endl;

  vector_t counts = posterior_atleast_one(dataset, present, absent);

  double m = 0;
  for(auto &x: counts)
    m += x;

  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "absent  = " << absent << endl
      << "m = " << m << endl;
  return(m);
};

HMM::posterior_t HMM::posterior_atleast_one(const Data::Seq &seq, bitmask_t present, bitmask_t absent) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq)"
      << "present = " << present << endl
      << "absent  = " << absent << endl;

  // TODO FIX ABSENT - done?

  double logp = log_likelihood_from_scale(compute_forward_scale(seq));

  double z;
  if(absent == 0) {
    SubHMM subhmm(*this, complementary_states_mask(present));
    double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(seq));

    z = 1 - exp(logp_wo_motif - logp);
    if(verbosity >= Verbosity::debug)
      cout << "seq = " << seq.definition << " " << seq.sequence << " logp = " << logp << " logp_wo_motif = " << logp_wo_motif << " z = " << z << endl;
  } else {
    SubHMM subhmm_wo_abs(*this, complementary_states_mask(absent));
    SubHMM subhmm_wo_abs_wo_motif(*this, complementary_states_mask(present | absent));
    if(verbosity >= Verbosity::debug) {
      cout << "Full     :" << *this << endl;
      cout << "Reduced 1:" << subhmm_wo_abs << endl;
      cout << "Reduced 2:" << subhmm_wo_abs_wo_motif << endl;
    }
    double logp_wo_abs = log_likelihood_from_scale(subhmm_wo_abs.compute_forward_scale(seq));
    double logp_wo_abs_wo_motif = log_likelihood_from_scale(subhmm_wo_abs_wo_motif.compute_forward_scale(seq));

    z = exp(logp_wo_abs - logp) - exp(logp_wo_abs_wo_motif - logp);
    if(verbosity >= Verbosity::debug) {
      stringstream s;
      s << "seq = " << seq.definition
        /*  << " " << seq.sequences[i].sequence */
        << " logp = " << logp
        << " logp_wo_abs = " << logp_wo_abs
        << " logp_wo_abs_wo_motif = " << logp_wo_abs_wo_motif
        << " z = " << z << endl;
      cout << s.str();
    }
  }

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq)" << endl
      << "present = " << present << endl
      << "absent  = " << absent << endl
      << "z = " << z << endl;
  posterior_t res = {logp, z};
  return(res);
};

double HMM::compute_score_all_motifs(const Data::Collection &collection, const Measures::Continuous::Measure &measure, bool weighting) const
{
  vector<size_t> all_motifs;
  for(size_t i = 0; i < groups.size(); i++)
    if(groups[i].kind == Group::Kind::Motif)
      all_motifs.push_back(i);
  return(compute_score(collection, measure, weighting, all_motifs, vector<size_t>()));
}

double HMM::compute_score(const Data::Collection &collection, const Measures::Continuous::Measure &measure, bool weighting, const vector<size_t> &present_motifs, const vector<size_t> &absent_motifs) const
{
  double score = 0;
  double W = 0;
  double w;
  switch(measure) {
    case Measure::Likelihood:
      score = log_likelihood(collection);
      break;
    case Measure::MutualInformation:
      for(auto &contrast: collection) {
        W += w = contrast.set_size;
        score += mutual_information(contrast, present_motifs, absent_motifs) * (weighting ? w : 1);
      }
      if(weighting)
        score /= W;
      break;
    case Measure::RankInformation:
      for(auto &contrast: collection)
        score = rank_information(contrast, present_motifs, absent_motifs);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      for(auto &contrast: collection)
        score += matthews_correlation_coefficient(contrast, present_motifs, absent_motifs);
      break;
    case Measure::DeltaFrequency:
      for(auto &contrast: collection)
        score += dips_sitescore(contrast, present_motifs, absent_motifs);
      break;
    case Measure::LogLikelihoodDifference:
      for(auto &contrast: collection)
        score += log_likelihood_difference(contrast, present_motifs, absent_motifs);
      break;
    case Measure::ClassificationPosterior:
      for(auto &contrast: collection)
        score += class_likelihood(contrast, present_motifs, absent_motifs, true);
      break;
    case Measure::ClassificationLikelihood:
      for(auto &contrast: collection)
        score += class_likelihood(contrast, present_motifs, absent_motifs, false);
      break;
    default:
      cout << "Score calculation for '" << measure2string(measure) << "' is not implemented." << endl;
      assert(0);
  }
  return(score);
}

double HMM::class_likelihood(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups, bool compute_posterior) const
{
  return(class_likelihood(contrast, make_mask(present_groups), make_mask(absent_groups), compute_posterior));
}

double HMM::class_likelihood(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent, bool compute_posterior) const
{
  double l = 0;
  for(auto &dataset: contrast)
    l += class_likelihood(dataset, present, absent, compute_posterior);
  if(verbosity >= Verbosity::debug)
    cout << "Data::Contrast l = " << l << endl;
  return(l);
}


double HMM::class_likelihood(const Data::Set &dataset, bitmask_t present, bitmask_t absent, bool compute_posterior) const
{
  // TODO FIX ABSENT
  // const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
  const double marginal_motif_prior  = 0.5;
  // TODO FIX ABSENT
  // const double class_cond = get_class_motif_prior(dataset.sha1, group_idx);
  const double class_cond = 0.2;
  const double log_class_prior = log(get_class_prior(dataset.sha1));

  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < dataset.set_size; i++) {
    posterior_t res = posterior_atleast_one(dataset.sequences[i], present, absent);
    double p = res.posterior;
    double x = log_class_prior + log(p * class_cond / marginal_motif_prior + (1-p) * (1-class_cond) / (1-marginal_motif_prior));
    if(not compute_posterior)
      x += res.log_likelihood;
    if(verbosity >= Verbosity::debug)
      cout << "Sequence " << dataset.sequences[i].definition << " p = " << p << " class log likelihood = " << x << " exp -> " << exp(x) << endl;
    l += x;
  }
  if(verbosity >= Verbosity::debug)
    cout << "Data::Set " << dataset.path << " l = " << l << endl;
  return(l);
}

double HMM::log_likelihood_difference(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  return(log_likelihood_difference(contrast, make_mask(present_groups), make_mask(absent_groups)));
}

double HMM::log_likelihood_difference(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  double d = 0;
  for(size_t sample_idx = 0; sample_idx < contrast.sets.size(); sample_idx++) {
    // TODO FIX ABSENT - done?
    bool signal = false;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if((bitmask_t(1 << group_idx) & present) != 0 and
          contrast.sets[sample_idx].motifs.find(groups[group_idx].name) != end(contrast.sets[sample_idx].motifs)) {
        signal = true;
        break;
      }
    d += (signal ? 1 : -1) * log_likelihood(contrast.sets[sample_idx]);
  }
  return(d);
}

double HMM::dips_sitescore(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &absent_groups) const
{
  return(dips_sitescore(contrast, make_mask(present_groups), make_mask(absent_groups)));
}

double HMM::dips_sitescore(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const
{
  vector_t posterior = posterior_atleast_one(contrast, present, absent);
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present, contrast, groups, false) + pseudo_count;
  size_t signal_size = m.true_positives + m.false_negatives;
  size_t control_size = m.false_positives + m.true_negatives;

  double t = m.true_positives / signal_size - m.false_positives / control_size;
  return(t);
}

double HMM::dips_tscore(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  return(dips_tscore(contrast, make_mask(present_groups)));
}

double HMM::dips_tscore(const Data::Contrast &contrast, bitmask_t present) const
{
  // TODO FIX ABSENT - done?
  vector_t posterior = expected_posterior(contrast, present);
  // TODO FIX ABSENT: adapt reduce - done?
  confusion_matrix m = reduce(posterior, present, contrast, groups, true) + pseudo_count;
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

