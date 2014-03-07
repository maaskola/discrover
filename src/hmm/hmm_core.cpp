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
 *       Filename:  hmm_core.cpp
 *
 *    Description:  Core routines of the HMM class 
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <boost/range/adaptors.hpp>
#include "../aux.hpp"
#include "hmm.hpp"

using namespace std;

double HMM::viterbi(const Data::Seq &s, StatePath &path) const
{
  // TODO: precompute log probabilities
  size_t L = s.isequence.size();

  vector_t v_current = scalar_vector(n_states, -numeric_limits<double>::infinity());
  vector_t v_previous = scalar_vector(n_states, -numeric_limits<double>::infinity());
  v_previous(start_state) = 0;
  boost::numeric::ublas::matrix<size_t> traceback(L,n_states);
  for(size_t i = 0; i < L; i++) {
    size_t symbol = s.isequence(i);
    if(symbol == empty_symbol)
      for(size_t k = start_state; k < n_states; k++) {
        double tmp = v_previous(k) + log(transition(k,start_state));
        if(tmp > v_current(start_state)) {
          v_current(start_state) = tmp;
          traceback(i,start_state) = k;
        }
      }
    else
      for(size_t l = start_state; l < n_states; l++) {
        double m = -numeric_limits<double>::infinity();
        for(size_t k = start_state; k < n_states; k++) {
          double tmp = v_previous(k) + log(transition(k,l));
          if(tmp > m) {
            m = tmp;
            traceback(i,l) = k;
          }
        }
        v_current(l) = log(emission(l, symbol)) + m;
      }
    v_previous = v_current;
    v_current = scalar_vector(n_states, -numeric_limits<double>::infinity());
  }

  double p = -numeric_limits<double>::infinity();
  size_t pi = 0;
  for(size_t k = start_state; k < n_states; k++) {
    double tmp = v_previous(k) + log(transition(k,start_state));
    if(tmp > p) {
      p = tmp;
      pi = k;
    }
  }
  path = StatePath(L);
  path(L-1) = pi;
  for(size_t i = L-1; i > 0; i--)
    pi = path(i-1) = traceback(i, pi);

  return(p);
};

vector_t HMM::compute_forward_scale(const Data::Seq &s) const
{
  size_t T = s.isequence.size();
  vector_t scale = zero_vector(T+2);
  vector_t prev = zero_vector(n_states);
  vector_t cur = zero_vector(n_states);

  prev(start_state) = 1;
  scale(0) = 1;
  for(size_t t = 0; t < T; t++) {
    size_t symbol = s.isequence(t);
    if(symbol == empty_symbol) {
      for(auto pre : pred[start_state])
        cur(start_state) += prev(pre) * transition(pre,start_state);
      scale(t+1) = cur(start_state);
      cur(start_state) = 1;
    } else {
      for(size_t i = 0; i < n_states; i++) {
        double emission_i_t = emission(i, symbol);
        if(emission_i_t > 0) {
          for(auto pre : pred[i])
            cur(i) += prev(pre) * transition(pre,i);
          cur(i) *= emission_i_t;
          scale(t+1) += cur(i);
        }
      }
      for(size_t i = 0; i < n_states; i++)
        cur(i) /= scale(t+1);
    }
    prev = cur;
    cur = zero_vector(n_states);
  }

  for(auto pre: pred[start_state])
    cur(start_state) += prev(pre) * transition(pre,start_state);
  scale(T+1) = cur(start_state);
  return(scale);
}

matrix_t HMM::compute_forward_scaled(const Data::Seq &s, vector_t &scale) const
{
  size_t T = s.isequence.size();
  matrix_t m = zero_matrix(T + 2, n_states);
  if(scale.size() != T + 2)
    scale = zero_vector(T+2);

  m(0,start_state) = 1;
  scale(0) = 1;
  for(size_t t = 0; t < T; t++) {
    size_t symbol = s.isequence(t);
    if(symbol == empty_symbol) {
      for(auto pre : pred[start_state])
        m(t+1,start_state) += m(t,pre) * transition(pre,start_state);
      scale(t+1) = m(t+1,start_state);
      m(t+1,start_state) = 1;
    } else {
      for(size_t i = 0; i < n_states; i++) {
        double emission_i_t = emission(i, symbol);
        if(emission_i_t > 0) {
          for(auto pre : pred[i])
            m(t+1,i) += m(t,pre) * transition(pre,i);
          m(t+1,i) *= emission_i_t;
          scale(t+1) += m(t+1,i);
        }
      }
      for(size_t i = 0; i < n_states; i++)
        m(t+1,i) /= scale(t+1);
    }
  }

  for(auto pre: pred[start_state])
    m(T+1,start_state) += m(T,pre) * transition(pre,start_state);
  scale(T+1) = m(T+1,start_state);
  m(T+1,start_state) = 1;

  if(verbosity >= Verbosity::debug)
    cout << "alpha = " << m << endl;
  return(m);
}

matrix_t HMM::compute_forward_prescaled(const Data::Seq &s, const vector_t &scale) const
{
  size_t T = s.isequence.size();
  matrix_t m = zero_matrix(T + 2, n_states);
  m(0,start_state) = 1.0/scale(0);
  for(size_t t = 0; t < T; t++) {
    size_t symbol = s.isequence(t);
    if(symbol == empty_symbol) {
      for(auto k: pred[start_state])
        m(t+1,start_state) += m(t,k) * transition(k,start_state);
      m(t+1,start_state) /= scale(t+1);
    } else {
      for(size_t i = 0; i < n_states; i++) {
        double emission_i_t = emission(i, symbol);
        if(emission_i_t > 0) {
          for(auto k: pred[i])
            m(t+1,i) += m(t,k) * transition(k,i);
          m(t+1,i) *= emission_i_t;
        }
        m(t+1,i) /= scale(t+1);
      }
    }
  }

  for(auto k: pred[start_state])
    m(T+1,start_state) += m(T,k) * transition(k,start_state);
  m(T+1,start_state) /= scale(T+1);
  return(m);
}

// Assuming that max_order == 0
matrix_t HMM::compute_backward_prescaled(const Data::Seq &s, const vector_t &scale) const
{
  size_t T = s.isequence.size();
  matrix_t m = zero_matrix(T + 2, n_states);
  m(T+1,start_state) = 1/scale(T+1);
  for(size_t i = 0; i < n_states; i++) {
    for(auto suc: succ[i])
      m(T,i) += m(T+1,suc) * transition(i,suc);
    m(T,i) /= scale(T);
  }

  for(int t = T-1; t >= 0; t--) {
    size_t symbol = s.isequence(t);
    if(symbol == empty_symbol)
      for(auto pre: pred[start_state])
        m(t,pre) = m(t+1,start_state) * transition(pre,start_state) / scale(t);
    else
      for(size_t i = 0; i < n_states; i++) { // TODO flip loops around
        for(auto suc: succ[i])
          m(t,i) += m(t+1,suc) * transition(i,suc) * emission(suc, symbol);
        m(t,i) /= scale(t);
      }
  }

  if(verbosity >= Verbosity::debug)
    cout << "beta = " << m << endl;
  return(m);
}

double HMM::likelihood_from_scale(const vector_t &scale) const
{
  double pf = 1;
  for(size_t i = 0; i < scale.size(); i++)
    pf *= scale(i);
  return(pf);
}

double HMM::log_likelihood_from_scale(const vector_t &scale) const
{
  double logpf = 0;
  for(size_t i = 0; i < scale.size(); i++)
    logpf += log(scale(i));
  return(logpf);
}

double HMM::expected_state_posterior(size_t k, const matrix_t &f, const matrix_t &b, const vector_t &scale) const
{
  double s = 0;
  for(size_t i = 0; i < f.size1(); i++)
    s += f(i,k) * b(i, k) * scale(i);
  return(s);
}

