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
#include "boost/multi_array.hpp"

using namespace std;

#define DO_PARALLEL 1

const bool verbose_conditional_mico_output = true;

HMM::pair_posterior_t &operator+=(HMM::pair_posterior_t &one, const HMM::pair_posterior_t &two)
{
  one.log_likelihood += two.log_likelihood;
  one.posterior_first += two.posterior_first;
  one.posterior_second += two.posterior_second;
  one.posterior_both += two.posterior_both;
  one.posterior_none += two.posterior_none;
  return(one);
}

ostream &operator<<(ostream &os, const HMM::pair_posterior_t &p)
{
  os << "logL = " << p.log_likelihood
    << " A = " << p.posterior_first
    << " B = " << p.posterior_second
    << " AB = " << p.posterior_both
    << " not(AB) = " << p.posterior_none;
  return(os);
}

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

double HMM::mutual_information(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing mutual information." << endl;
  bitmask_t present = make_mask(present_groups);
  return(mutual_information(contrast, present));
}

double HMM::conditional_mutual_information(const Data::Contrast &contrast, const vector<size_t> &present_groups, const vector<size_t> &previous_groups, double residual_ratio_threshold) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing residual mutual information." << endl;
  bitmask_t present = make_mask(present_groups);
  bitmask_t previous = make_mask(previous_groups);
  return(conditional_mutual_information(contrast, present, previous, residual_ratio_threshold));
}

double HMM::rank_information(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing mutual information." << endl;
  bitmask_t present = make_mask(present_groups);
  return(rank_information(contrast, present));
}

double HMM::matthews_correlation_coefficient(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  if(verbosity >= Verbosity::debug)
    cout << "Computing sum of Matthew's correlation coefficient." << endl;
  bitmask_t present = make_mask(present_groups);
  return(matthews_correlation_coefficient(contrast, present));
}

confusion_matrix HMM::reduce(const vector_t &v, HMM::bitmask_t present, const Data::Contrast &contrast, bool word_stats) const
{
  confusion_matrix m = {0, 0, 0, 0};
  for(size_t sample_idx = 0; sample_idx < v.size(); sample_idx++)
    if(is_present(contrast.sets[sample_idx], present)) {
      m.true_positives += v[sample_idx];
      m.false_negatives += (word_stats ? contrast.sets[sample_idx].seq_size : contrast.sets[sample_idx].set_size) - v[sample_idx];
    } else {
      m.false_positives += v[sample_idx];
      m.true_negatives += (word_stats ? contrast.sets[sample_idx].seq_size : contrast.sets[sample_idx].set_size) - v[sample_idx];
    }
  return(m);
}

double HMM::mutual_information(const Data::Contrast &contrast, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information(Data::Contrast)" << endl;
  vector_t posterior = posterior_atleast_one(contrast, present);
  matrix_t m(posterior.size(), 2);
  for(size_t i = 0; i < posterior.size(); i++) {
    m(i,0) = posterior(i);
    m(i,1) = contrast.sets[i].set_size - posterior(i);
  }
  m = m + pseudo_count;
  if(verbosity >= Verbosity::debug)
    cout << "HMM::mutual_information(Data::Contrast)" << endl
      << "present = " << present << endl
      << "counts  = " << m << endl;
  double mi = calc_mutual_information(m, 0, true, false, false);
//  if(not check_enrichment(contrast, m, group_idx))
//    mi = -mi;
  return(mi);
}

/** The conditional mutual information
 * This is the reduction in the uncertainty of X due to the knowledge of Y when Z is given:
 * I( X; Y | Z ) = H( X | Z ) - H( X | Y, Z )
 * Symmetrically, it is also the reduction in the uncertainty of Y due to the knowledge of X when Z is given:
 * I( X; Y | Z ) = H( Y | Z ) - H( Y | X, Z )
 *
 * Here X is the first motif, Y is the contrast, and Z is the second (previous) motif(s)
 *
 * See Cover & Thomas 2006 equations (2.60) and (2.61)
 */
double calc_conditional_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity)
{
  const size_t X = 2;
  const size_t Y = contrast.sets.size();
  const size_t Z = 2;

  typedef boost::multi_array<double, 1> array1_t;
  typedef array1_t::index index;

  typedef boost::multi_array<double, 2> array2_t;
  typedef array2_t::index index;

  typedef boost::multi_array<double, 3> array_t;
  typedef array_t::index index;
  //
  // the joint probability of X, Y, and Z
  array_t p(boost::extents[X][Y][Z]);

 // fill joint probability table with absolute counts
  for(size_t i = 0; i < contrast.sets.size(); i++) {
    p[0][i][0] = pair_posteriors[i].posterior_both;
    p[0][i][1] = pair_posteriors[i].posterior_first - pair_posteriors[i].posterior_both;
    p[1][i][0] = pair_posteriors[i].posterior_second - pair_posteriors[i].posterior_both;
    p[1][i][1] = pair_posteriors[i].posterior_none;
  }

  // add pseudo-count
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        p[x][y][z] += ps;

  double marginal = 0;
  // sum over entries to compute marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        marginal += p[x][y][z];

  // normalize joint probability by dividing through marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        p[x][y][z] /= marginal;


  // the joint marginal probability of X and Z
  array2_t pxz(boost::extents[X][Z]);
  // the joint marginal probability of Y and Z
  array2_t pyz(boost::extents[Y][Z]);

  // compute marginal distribution of X and Z by summing over Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pxz[x][z] += p[x][y][z];

  // compute marginal distribution of Y and Z by summing over X
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pyz[y][z] += p[x][y][z];


  // the marginal probability of Z
  array1_t pz(boost::extents[Z]);

  // compute marginal distribution of Z by summing over X and Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pz[z] += p[x][y][z];

  double mi = 0;
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        mi += p[x][y][z] * log( p[x][y][z] / pxz[x][z] / pyz[y][z] * pz[z]);

  mi /= log(2.0);

  mi = max<double>(mi, 0);

  if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info)) {
    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information z = " << z << " p(x,y,z) =";
      for(size_t x = 0; x < X; x++)
        for(size_t y = 0; y < Y; y++)
          cout << " " << p[x][y][z];
      cout << endl;
    }

    cout << "conditional_mutual_information p(x,z) =";
    for(size_t x = 0; x < X; x++)
      for(size_t z = 0; z < Z; z++)
        cout << " " << pxz[x][z];
    cout << endl;

    cout << "conditional_mutual_information p(y,z) =";
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        cout << " " << pyz[y][z];
    cout << endl;

    cout << "conditional_mutual_information pz =";
    for(size_t z = 0; z < Z; z++)
      cout << " " << pz[z];
    cout << endl;

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information z = " << z << " p(x,y|z) =";
      for(size_t x = 0; x < X; x++)
        for(size_t y = 0; y < Y; y++)
          cout << " " << p[x][y][z] / pz[z];
      cout << endl;
    }

    cout << "conditional_mutual_information p(x|z) =";
    for(size_t x = 0; x < X; x++)
      for(size_t z = 0; z < Z; z++)
        cout << " " << pxz[x][z] / pz[z];
    cout << endl;

    cout << "conditional_mutual_information p(y|z) =";
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        cout << " " << pyz[y][z] / pz[z];
    cout << endl;
  }

  return(mi);
}

/** The pair mutual information
 * This is the mutual information of the two motifs
 */
double pair_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity)
{
  const size_t X = 2;
  const size_t Y = 2;

  typedef boost::multi_array<double, 1> array1_t;
  typedef array1_t::index index;

  typedef boost::multi_array<double, 2> array2_t;
  typedef array2_t::index index;

  // the joint probability of X and Y
  array2_t p(boost::extents[X][Y]);

 // fill joint probability table with absolute counts
  for(size_t i = 0; i < contrast.sets.size(); i++) {
    p[0][0] += pair_posteriors[i].posterior_both;
    p[0][1] += pair_posteriors[i].posterior_first - pair_posteriors[i].posterior_both;
    p[1][0] += pair_posteriors[i].posterior_second - pair_posteriors[i].posterior_both;
    p[1][1] += pair_posteriors[i].posterior_none;
  }

  // add pseudo-count
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      p[x][y] += ps;

  double marginal = 0;
  // sum over entries to compute marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      marginal += p[x][y];

  // normalize joint probability by dividing through marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      p[x][y] /= marginal;


  // the marginal probability of X
  array1_t px(boost::extents[X]);
  // the marginal probability of Y
  array1_t py(boost::extents[Y]);

  // compute marginal distribution of X by summing over Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
        px[x] += p[x][y];

  // compute marginal distribution of Y by summing over X
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
        py[y] += p[x][y];

  if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info)) {
    cout << "pair_mutual_information joint =";
    for(size_t x = 0; x < X; x++)
      for(size_t y = 0; y < Y; y++)
        cout << " " << p[x][y];
    cout << endl;
    cout << "pair_mutual_information px =";
    for(size_t x = 0; x < X; x++)
      cout << " " << px[x];
    cout << endl;
    cout << "pair_mutual_information py =";
    for(size_t y = 0; y < Y; y++)
      cout << " " << py[y];
    cout << endl;
  }

  double mi = 0;
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      mi += p[x][y] * log( p[x][y] / px[x] / py[y] );

  mi /= log(2.0);

  mi = max<double>(mi, 0);

  return(mi);
}

pair<double,double> HMM::conditional_and_motif_pair_mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::conditional_and_motif_pair_mutual_information(Data::Contrast)" << endl;
  auto pair_posteriors = pair_posterior_atleast_one(contrast, present, previous);
  double conditional_mi = calc_conditional_mutual_information(contrast, pair_posteriors, pseudo_count, verbosity);
  double pair_mi = pair_mutual_information(contrast, pair_posteriors, pseudo_count, verbosity);
  if(verbosity >= Verbosity::verbose)
    cout << "HMM::conditional_and_motif_pair_mutual_information(Data::Contrast)" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl
      << "condMI   = " << conditional_mi << " pairMI = " << pair_mi << " ratio = " << (conditional_mi / pair_mi) << endl;
  return(make_pair(conditional_mi, pair_mi));
}


double HMM::conditional_mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::conditional_mutual_information(Data::Contrast)" << endl;
  auto pair_posteriors = pair_posterior_atleast_one(contrast, present, previous);
  double conditional_mi = calc_conditional_mutual_information(contrast, pair_posteriors, pseudo_count, verbosity);
  if(verbosity >= Verbosity::verbose)
    cout << "HMM::conditional_mutual_information(Data::Contrast)" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl
      << "condMI   = " << conditional_mi << endl;
  return(conditional_mi);
}

double HMM::motif_pair_mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::motif_pair_mutual_information(Data::Contrast)" << endl;
  auto pair_posteriors = pair_posterior_atleast_one(contrast, present, previous);
  double pair_mi = pair_mutual_information(contrast, pair_posteriors, pseudo_count, verbosity);
  if(verbosity >= Verbosity::verbose)
    cout << "HMM::conditional_mutual_information(Data::Contrast)" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl
      << "pairMI = " << pair_mi << endl;
  return(pair_mi);
}

double HMM::rank_information(const Data::Contrast &contrast, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Contrast)" << endl;
  double ri = 0;
  for(auto &dataset: contrast)
    ri += rank_information(dataset, present);
  return(ri);
}


double HMM::rank_information(const Data::Set &dataset, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::rank_information(Data::Set)" << endl;
  vector_t posterior = posterior_atleast_one(dataset, present);
  return(calc_rank_information(posterior, pseudo_count));
}

double HMM::matthews_correlation_coefficient(const Data::Contrast &contrast, bitmask_t present) const
{
  // TODO find out if there's a proper generalization of the MCC to multiple experiments
  vector_t posterior = posterior_atleast_one(contrast, present);
  confusion_matrix m = reduce(posterior, present, contrast, false) + pseudo_count;
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

vector_t HMM::posterior_atleast_one(const Data::Contrast &contrast, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Contrast)" << endl
      << "present =" << present << endl;
  vector_t v(contrast.sets.size());
  for(size_t i = 0; i < contrast.sets.size(); i++)
    v[i] = sum_posterior_atleast_one(contrast.sets[i], present);
  return(v);
}

HMM::pair_posteriors_t HMM::pair_posterior_atleast_one(const Data::Contrast &contrast, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::pair_posterior_atleast_one(Data::Contrast)" << endl
      << "present  =" << present << endl
      << "previous = " << previous << endl;

  pair_posteriors_t v(contrast.sets.size());
  for(size_t i = 0; i < contrast.sets.size(); i++)
    v[i] = sum_pair_posterior_atleast_one(contrast.sets[i], present, previous);

  if(verbosity >= Verbosity::debug)
    // TODO put results into debug output
    cout << "HMM::pair_posterior_atleast_one(Data::Contrast)" << endl
      << "present  =" << present << endl
      << "previous = " << previous << endl;
  return(v);
}


vector_t HMM::posterior_atleast_one(const Data::Set &dataset, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl;

  if(verbosity >= Verbosity::debug) {
    cout << "complementary_states_mask(present)) =";
    for(auto &state: complementary_states_mask(present))
      cout << " " << state;
    cout << endl;
  }

  vector_t vec(dataset.set_size);
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

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "vec = " << vec << endl;
  return(vec);
};


HMM::pair_posteriors_t HMM::pair_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::pair_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl;

  if(verbosity >= Verbosity::debug) {
    cout << "complementary_states_mask(present)) =";
    for(auto &state: complementary_states_mask(present))
      cout << " " << state;
    cout << endl;
  }

  pair_posteriors_t vec(dataset.set_size);
  SubHMM subhmm_one(*this, complementary_states_mask(present));
  SubHMM subhmm_two(*this, complementary_states_mask(previous));
  SubHMM subhmm_both(*this, complementary_states_mask(present | previous));
#pragma omp parallel for schedule(static) if(DO_PARALLEL)
  for(size_t i = 0; i < dataset.set_size; i++) {
    double logp = log_likelihood_from_scale(compute_forward_scale(dataset.sequences[i]));
    double logp_wo_one = log_likelihood_from_scale(subhmm_one.compute_forward_scale(dataset.sequences[i]));
    double logp_wo_two = log_likelihood_from_scale(subhmm_two.compute_forward_scale(dataset.sequences[i]));
    double logp_wo_either = log_likelihood_from_scale(subhmm_both.compute_forward_scale(dataset.sequences[i]));
    double z_one = 1 - exp(logp_wo_one - logp);
    double z_two = 1 - exp(logp_wo_two - logp);
    double z_either = 1 - exp(logp_wo_either - logp);
    double z_both =  z_one + z_two - z_either;
    pair_posterior_t p = {logp, z_one, z_two, z_both, 1 - z_either};
    if(verbosity >= Verbosity::debug)
      cout << "seq = " << dataset.sequences[i].definition << p << endl;
    vec[i] = p;
  }

  if(verbosity >= Verbosity::debug)
    cout << "HMM::pair_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl;
  return(vec);
}

double HMM::sum_posterior_atleast_one(const Data::Set &dataset, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl;

  vector_t counts = posterior_atleast_one(dataset, present);

  double m = 0;
  for(auto &x: counts)
    m += x;

  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present = " << present << endl
      << "m = " << m << endl;
  return(m);
};

HMM::pair_posterior_t HMM::sum_pair_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t previous) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::sum_pair_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl;

  pair_posteriors_t pair_counts = pair_posterior_atleast_one(dataset, present, previous);

  pair_posterior_t summed_pair_counts = {0, 0, 0, 0, 0};
  for(auto &x: pair_counts)
    summed_pair_counts += x;

  if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info))
    cout << "HMM::sum_pair_posterior_atleast_one(Data::Set = " << dataset.path << ")" << endl
      << "present  = " << present << endl
      << "previous = " << previous << endl
      << "counts: " << summed_pair_counts << endl;
  return(summed_pair_counts);
}

HMM::posterior_t HMM::posterior_atleast_one(const Data::Seq &seq, bitmask_t present) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq)"
      << "present = " << present << endl;

  double logp = log_likelihood_from_scale(compute_forward_scale(seq));

  double z;
  SubHMM subhmm(*this, complementary_states_mask(present));
  double logp_wo_motif = log_likelihood_from_scale(subhmm.compute_forward_scale(seq));

  z = 1 - exp(logp_wo_motif - logp);
  if(verbosity >= Verbosity::debug)
    cout << "seq = " << seq.definition << " " << seq.sequence << " logp = " << logp << " logp_wo_motif = " << logp_wo_motif << " z = " << z << endl;

  if(verbosity >= Verbosity::debug)
    cout << "HMM::posterior_atleast_one(Data::Seq)" << endl
      << "present = " << present << endl
      << "z = " << z << endl;
  posterior_t res = {logp, z};
  return(res);
};

double HMM::compute_score_all_motifs(const Data::Collection &collection, const Measures::Continuous::Measure &measure, const Options::HMM &options) const
{
  vector<size_t> all_motifs;
  for(size_t i = 0; i < groups.size(); i++)
    if(groups[i].kind == Group::Kind::Motif)
      all_motifs.push_back(i);
  return(compute_score(collection, measure, options, all_motifs, vector<size_t>()));
}

double HMM::compute_score(const Data::Collection &collection, const Measures::Continuous::Measure &measure, const Options::HMM &options, const vector<size_t> &present_motifs, const vector<size_t> &previous_motifs) const
{
  double score = 0;
  double W = 0;
  double w;
  bitmask_t present = make_mask(present_motifs);
  bitmask_t previous = make_mask(previous_motifs);
  switch(measure) {
    case Measure::Likelihood:
      score = log_likelihood(collection);
      break;
    case Measure::MutualInformation:
      for(auto &contrast: collection) {
        W += w = contrast.set_size;
        score += mutual_information(contrast, present) * (options.weighting ? w : 1);
      }
      if(options.weighting)
        score /= W;
      break;
    case Measure::ConditionalPairMutualInformationRatio:
    case Measure::ThresholdedConditionalMutualInformation:
      {
        double conditional_mi = 0;
        double pair_mi = 0;
        for(auto &contrast: collection) {
          W += w = contrast.set_size;
          auto current = conditional_and_motif_pair_mutual_information(contrast, present, previous);
          conditional_mi += current.first * (options.weighting ? w : 1);
          pair_mi += current.second * (options.weighting ? w : 1);
        }
        if(options.weighting) {
          conditional_mi /= W;
          pair_mi /= W;
        }
        double ratio = conditional_mi / pair_mi;
        if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info))
          cout << "condMI   = " << conditional_mi << " pairMI = " << pair_mi << " ratio = " << ratio << endl;
        if(measure == Measure::ConditionalPairMutualInformationRatio)
          return(ratio);
        score = conditional_mi;
        if(ratio < options.multi_motif.residual_ratio)
          score = 0;
      }
      break;
    case Measure::ConditionalMutualInformation:
      for(auto &contrast: collection) {
        W += w = contrast.set_size;
        score += conditional_mutual_information(contrast, present, previous) * (options.weighting ? w : 1);
      }
      if(options.weighting)
        score /= W;
      break;
    case Measure::PairMutualInformation:
      for(auto &contrast: collection) {
        W += w = contrast.set_size;
        score += motif_pair_mutual_information(contrast, present, previous) * (options.weighting ? w : 1);
      }
      if(options.weighting)
        score /= W;
      break;
    case Measure::RankInformation:
      for(auto &contrast: collection)
        score = rank_information(contrast, present);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      for(auto &contrast: collection)
        score += matthews_correlation_coefficient(contrast, present);
      break;
    case Measure::DeltaFrequency:
      for(auto &contrast: collection)
        score += dips_sitescore(contrast, present);
      break;
    case Measure::LogLikelihoodDifference:
      for(auto &contrast: collection)
        score += log_likelihood_difference(contrast, present);
      break;
    case Measure::ClassificationPosterior:
      for(auto &contrast: collection)
        score += class_likelihood(contrast, present, true);
      break;
    case Measure::ClassificationLikelihood:
      for(auto &contrast: collection)
        score += class_likelihood(contrast, present, false);
      break;
    default:
      cout << "Score calculation for '" << measure2string(measure) << "' is not implemented." << endl;
      assert(0);
  }
  return(score);
}

double HMM::class_likelihood(const Data::Contrast &contrast, const vector<size_t> &present_groups, bool compute_posterior) const
{
  return(class_likelihood(contrast, make_mask(present_groups), compute_posterior));
}

double HMM::class_likelihood(const Data::Contrast &contrast, bitmask_t present, bool compute_posterior) const
{
  double l = 0;
  for(auto &dataset: contrast)
    l += class_likelihood(dataset, present, compute_posterior);
  if(verbosity >= Verbosity::debug)
    cout << "Data::Contrast l = " << l << endl;
  return(l);
}


double HMM::class_likelihood(const Data::Set &dataset, bitmask_t present, bool compute_posterior) const
{
  // TODO FIX BITMASK MMIE
  // const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
  const double marginal_motif_prior  = 0.5;
  // TODO FIX BITMASK MMIE
  // const double class_cond = get_class_motif_prior(dataset.sha1, group_idx);
  const double class_cond = 0.2;
  const double log_class_prior = log(get_class_prior(dataset.sha1));

  double l = 0;
#pragma omp parallel for schedule(static) reduction(+:l) if(DO_PARALLEL)
  for(size_t i = 0; i < dataset.set_size; i++) {
    posterior_t res = posterior_atleast_one(dataset.sequences[i], present);
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

double HMM::log_likelihood_difference(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  return(log_likelihood_difference(contrast, make_mask(present_groups)));
}

bool HMM::is_present(const Data::Set &dataset, bitmask_t present) const
{
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if((bitmask_t(1 << group_idx) & present) != 0 and
        dataset.motifs.find(groups[group_idx].name) != end(dataset.motifs))
      return(true);
  return(false);
}

double HMM::log_likelihood_difference(const Data::Contrast &contrast, bitmask_t present) const
{
  double d = 0;
  for(size_t sample_idx = 0; sample_idx < contrast.sets.size(); sample_idx++) {
    bool signal = is_present(contrast.sets[sample_idx], present);
    d += (signal ? 1 : -1) * log_likelihood(contrast.sets[sample_idx]);
  }
  return(d);
}

double HMM::dips_sitescore(const Data::Contrast &contrast, const vector<size_t> &present_groups) const
{
  return(dips_sitescore(contrast, make_mask(present_groups)));
}

double HMM::dips_sitescore(const Data::Contrast &contrast, bitmask_t present) const
{
  vector_t posterior = posterior_atleast_one(contrast, present);
  confusion_matrix m = reduce(posterior, present, contrast, false) + pseudo_count;
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
  vector_t posterior = expected_posterior(contrast, present);
  confusion_matrix m = reduce(posterior, present, contrast, true) + pseudo_count;
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

