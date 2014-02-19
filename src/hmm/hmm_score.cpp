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

#define DO_PARALLEL 1

double HMM::mutual_information(const Data::Series &data, const std::vector<size_t> &groups) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Computing feature-wise mutual information." << std::endl;
  double mi = 0;
  for(auto &group_idx: groups)
    mi += mutual_information(data, group_idx);
  return(mi);
}

double HMM::rank_information(const Data::Series &data, const std::vector<size_t> &groups) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Computing feature-wise mutual information." << std::endl;
  double mi = 0;
  for(auto &group_idx: groups)
    mi += rank_information(data, group_idx);
  return(mi);
}

double HMM::matthews_correlation_coefficient(const Data::Series &data, const std::vector<size_t> &groups) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Computing sum of feature-wise Matthew's correlation coefficient." << std::endl;
  double mcc = 0;
  for(auto &group_idx: groups)
    mcc += matthews_correlation_coefficient(data, group_idx);
  return(mcc);
}

confusion_matrix reduce(const vector_t &v, size_t group_idx, const Data::Series &data, const std::vector<Group> &groups, bool word_stats) {
  confusion_matrix m = {0, 0, 0, 0};
  for(size_t sample_idx = 0; sample_idx < v.size(); sample_idx++) {
    if(std::find(data.sets[sample_idx].motifs.begin(), data.sets[sample_idx].motifs.end(), groups[group_idx].name) != data.sets[sample_idx].motifs.end()) {
      m.true_positives += v[sample_idx];
      m.false_negatives += (word_stats ? data.sets[sample_idx].seq_size : data.sets[sample_idx].set_size) - v[sample_idx];
    } else {
      m.false_positives += v[sample_idx];
      m.true_negatives += (word_stats ? data.sets[sample_idx].seq_size : data.sets[sample_idx].set_size) - v[sample_idx];
    }
  }
  return(m);
}

double HMM::mutual_information(const Data::Series &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::mutual_information(Data::Series, Feature)" << std::endl;
  vector_t posterior = posterior_atleast_one(data, group_idx);
  matrix_t m(posterior.size(), 2);
  for(size_t i = 0; i < posterior.size(); i++) {
    m(i,0) = posterior(i);
    m(i,1) = data.sets[i].set_size - posterior(i);
  }
  m = m + pseudo_count;
  double mi = calc_mutual_information(m, 0, true, false, false);
//  if(not check_enrichment(data, m, group_idx))
//    mi = -mi;
  return(mi);
}

double HMM::rank_information(const Data::Series &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::rank_information(Data::Series, Feature)" << std::endl;
  double ri = 0;
  for(auto &set: data)
    ri += rank_information(set, group_idx);
  return(ri);
}


double HMM::rank_information(const Data::Set &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::rank_information(Data::Set, Feature)" << std::endl;
  vector_t posterior = posterior_atleast_one(data, group_idx);
  return(calc_rank_information(posterior, pseudo_count));
}

double HMM::matthews_correlation_coefficient(const Data::Series &data, size_t group_idx) const
{
  // TODO find out if there's a proper generalization of the MCC to multiple experiments
  vector_t posterior = posterior_atleast_one(data, group_idx);
  confusion_matrix m = reduce(posterior, group_idx, data, groups, false) + pseudo_count;
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

vector_t HMM::posterior_atleast_one(const Data::Series &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::posterior_atleast_one(Data::Series, group_idx=" << group_idx << ")" << std::endl;
  vector_t v(data.sets.size());
  for(size_t i = 0; i < data.sets.size(); i++)
    v[i] = sum_posterior_atleast_one(data.sets[i], group_idx);
  return(v);
}

vector_t HMM::posterior_atleast_one(const Data::Set &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::posterior_atleast_one(Data::Set = " << data.path << ", group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << std::endl;

  SubHMM subhmm(*this, complementary_states(group_idx));

  vector_t vec(data.set_size);
  // std::cout << subhmm << std::endl;
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++) {
    double logp = log_likelihood_from_scale(compute_forward_scale(data.sequences[i]));
    double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(data.sequences[i]));

    double z = 1 - exp(logp_wo_motif - logp);
    if(verbosity >= Verbosity::debug)
      std::cout << "seq = " << data.sequences[i].definition << " " << data.sequences[i].sequence << " logp = " << logp << " logp_wo_motif = " << logp_wo_motif << " z = " << z << std::endl;
    vec[i] = z;
  }

  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::posterior_atleast_one(Data::Set = " << data.path << ", group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << " vec = " << vec << std::endl;
  return(vec);
};



double HMM::sum_posterior_atleast_one(const Data::Set &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::sum_posterior_atleast_one(Data::Set = " << data.path << ", group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << std::endl;

  vector_t counts = posterior_atleast_one(data, group_idx);

  double m = 0;
  for(auto &x: counts)
    m += x;

  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::sum_posterior_atleast_one(Data::Set = " << data.path << ", group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << " m = " << m << std::endl;
  return(m);
};

HMM::posterior_t HMM::posterior_atleast_one(const Data::Seq &data, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::posterior_atleast_one(Data::Seq, group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << std::endl;

  SubHMM subhmm(*this, complementary_states(group_idx));

  // std::cout << subhmm << std::endl;
  double logp = log_likelihood_from_scale(compute_forward_scale(data));
  double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(data));

  double z = 1 - exp(logp_wo_motif - logp);
  if(verbosity >= Verbosity::debug)
    std::cout << "seq = " << data.definition << " " << data.sequence << " logp = " << logp << " logp_wo_motif = " << logp_wo_motif << " z = " << z << std::endl;

  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::posterior_atleast_one(Data::Seq, group_idx = " << group_idx << " - '" << groups[group_idx].name << "')" << " z = " << z << std::endl;
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
    if(std::find(path.begin(), path.end(), groups[group_idx].states[0]) != path.end())
      m++;
  }
  return(m);
};

double HMM::compute_score(const Data::Collection &data, const Training::Task &task, bool weighting) const
{
  std::vector<size_t> all_motifs;
  for(size_t i = 0; i < groups.size(); i++)
    if(groups[i].kind == Group::Kind::Motif)
      all_motifs.push_back(i);
  double score = 0;
  double W = 0;
  double w;
  switch(task.measure) {
    case Measure::Likelihood:
      score = log_likelihood(data);
      break;
    case Measure::MutualInformation:
      for(auto &series: data) {
        W += w = series.set_size;
        score += mutual_information(series, all_motifs) * (weighting ? w : 1);
      }
      if(weighting)
        score /= W;
      break;
    case Measure::RankInformation:
      for(auto &series: data)
        score = rank_information(series, all_motifs);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      for(auto &series: data)
        score += matthews_correlation_coefficient(series, all_motifs);
      break;
    case Measure::DeltaFrequency:
      for(auto &series: data)
        score += dips_sitescore(series, all_motifs);
      break;
    case Measure::LogLikelihoodDifference:
      for(auto &series: data)
        score += log_likelihood_difference(series, all_motifs);
      break;
    case Measure::ClassificationPosterior:
      for(auto &series: data)
        score += class_likelihood(series, all_motifs, true);
      break;
    case Measure::ClassificationLikelihood:
      for(auto &series: data)
        score += class_likelihood(series, all_motifs, false);
      break;
    default:
      std::cout << "Score calculation for '" << measure2string(task.measure) << "' is not implemented." << std::endl;
      assert(0);
  }
  return(score);
}

double HMM::class_likelihood(const Data::Series &data, const std::vector<size_t> &groups, bool compute_posterior) const
{
  double l = 0;
  for(auto group_idx: groups)
    l += class_likelihood(data, group_idx, compute_posterior);
  return(l);
}

double HMM::class_likelihood(const Data::Series &data, size_t group_idx, bool compute_posterior) const
{
  double l = 0;
  for(auto &set: data)
    l += class_likelihood(set, group_idx, compute_posterior);
  if(verbosity >= Verbosity::debug)
    std::cout << "Data::Series l = " << l << std::endl;
  return(l);
}


double HMM::class_likelihood(const Data::Set &data, size_t group_idx, bool compute_posterior) const
{
  const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
  const double class_cond = get_class_motif_prior(data.sha1, group_idx);
  const double log_class_prior = log(get_class_prior(data.sha1));

  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < data.set_size; i++) {
    posterior_t res = posterior_atleast_one(data.sequences[i], group_idx);
    double p = res.posterior;
    double x = log_class_prior + log(p * class_cond / marginal_motif_prior + (1-p) * (1-class_cond) / (1-marginal_motif_prior));
    if(not compute_posterior)
      x += res.log_likelihood;
    if(verbosity >= Verbosity::debug)
      std::cout << "Sequence " << data.sequences[i].definition << " p = " << p << " class log likelihood = " << x << " exp -> " << exp(x) << std::endl;
    l += x;
  }
  if(verbosity >= Verbosity::debug)
    std::cout << "Data::Set " << data.path << " l = " << l << std::endl;
  return(l);
}

double HMM::log_likelihood_difference(const Data::Series &data, const std::vector<size_t> &groups) const
{
  double s = 0;
  for(auto group_idx: groups)
    s += log_likelihood_difference(data, group_idx);
  return(s);
}

double HMM::log_likelihood_difference(const Data::Series &data, size_t group_idx) const
{
  double d = 0;
  for(size_t sample_idx = 0; sample_idx < data.sets.size(); sample_idx++) {
    bool signal = false;
    if(data.sets[sample_idx].motifs.find(groups[group_idx].name) != data.sets[sample_idx].motifs.end())
      signal = true;
    d += (signal ? 1 : -1) * log_likelihood(data.sets[sample_idx]);
  }
  return(d);
}

double HMM::dips_sitescore(const Data::Series &data, const std::vector<size_t> &groups) const
{
  double t = 0;
  for(auto &group_idx: groups)
    t += dips_sitescore(data, group_idx);
  return(t);
}

double HMM::dips_sitescore(const Data::Series &data, size_t group_idx) const
{
  vector_t posterior = posterior_atleast_one(data, group_idx);
  confusion_matrix m = reduce(posterior, group_idx, data, groups, false) + pseudo_count;
  size_t signal_size = m.true_positives + m.false_negatives;
  size_t control_size = m.false_positives + m.true_negatives;

  double t = m.true_positives / signal_size - m.false_positives / control_size;
  return(t);
}

double HMM::dips_tscore(const Data::Series &data, const std::vector<size_t> &groups) const
{
  double t = 0;
  for(auto &group_idx: groups)
    t += dips_tscore(data, group_idx);
  return(t);
}

double HMM::dips_tscore(const Data::Series &data, size_t group_idx) const
{
  vector_t posterior = expected_posterior(data, group_idx);
  confusion_matrix m = reduce(posterior, group_idx, data, groups, true) + pseudo_count;
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

