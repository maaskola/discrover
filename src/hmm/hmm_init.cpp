/* =====================================================================================
 * Copyright (c) 2012, Jonas Maaskola
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
 *       Filename:  hmm_init.cpp
 *
 *    Description:  Auxiliary routines to initialize HMMs.
 *
 *        Created:  05/11/2012 03:59:10 AM
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include "hmm.hpp"

void HMM::finalize_initialization()
{
  initialize_ranges();
  initialize_pred_succ();
  check_consistency();
}

void HMM::initialize_ranges()
{
  all_range = Training::Range();
  constitutive_range = Training::Range();
  emitting_range = Training::Range();
  for(size_t i = start_state; i < n_states; i++) {
    all_range.push_back(i);
    if(group_ids[i] == 0) // TODO this could use more information that is available in data structure Group::Kind.
      constitutive_range.push_back(i);
    if(i != start_state)
      emitting_range.push_back(i);
  }
  if(verbosity >= Verbosity::debug) {
    std::cout << "all_range =";
    for(auto x: all_range)
      std::cout << " " << x;
    std::cout << std::endl;
    std::cout << "constitutive_range =";
    for(auto x: constitutive_range)
      std::cout << " " << x;
    std::cout << std::endl;
    std::cout << "emitting_range =";
    for(auto x: emitting_range)
      std::cout << " " << x;
    std::cout << std::endl;
  }
}

void HMM::initialize_pred_succ()
{
  pred = std::vector<std::list<size_t>>();
  succ = std::vector<std::list<size_t>>();
  for(size_t i = 0; i < n_states; i++) {
    pred.push_back(std::list<size_t>());
    succ.push_back(std::list<size_t>());
  }
  for(size_t i = 0; i < n_states; i++)
    for(size_t j = 0; j < n_states; j++)
      if(transition(i,j) > 0) {
        pred[j].push_back(i);
        succ[i].push_back(j);
      }
};

/** Initialize the emission matrix */
void HMM::initialize_emissions(size_t bg_order)
{
  // Start state does not emit
  for(size_t i = 0; i < n_emissions; i++)
    emission(start_state,i) = 0;

  // All other states are uniform over all emissions
  for(size_t i = bg_state; i < n_states; i++)
    for(size_t j = 0; j < n_emissions; j++)
      emission(i,j) = 1.0 / n_emissions;
  // TODO set up the emissions for higher order emission states
};

/** Initialize the transition matrix to zero */
void HMM::initialize_transitions()
{
  transition = zero_matrix(n_states, n_states);
}

void HMM::initialize_bg_transitions()
{
  transition(start_state, bg_state) = 1.0;
  transition(bg_state, start_state) = 1.0;
  transition(bg_state, bg_state) = 1.0;
};


/** Initialize the transition probabilities to and from the motif chain such
  * that lambda of the sequences have a motif occurrence, where the expected
  * sequence length is l.
  */
void HMM::initialize_transitions_to_and_from_chain(size_t w, double l, double lambda, size_t first, size_t last, size_t pad_left, size_t pad_right)
{
  const double x = lambda / (l - w + 1);
  const double y = lambda / (l - w + 1) * (l - w) / (l - w * lambda);
  const double z = (1 - lambda / (l - w + 1)) / (l - w * lambda);

  if(verbosity >= Verbosity::debug)
    std::cout << "l = " << l << std::endl
      << "w = " << w << std::endl
      << "x = " << x << std::endl
      << "y = " << y << std::endl
      << "z = " << z << std::endl;

  // initialize transitions from the start state
  transition(start_state,start_state) = 0; // sequences are at least one position long
  transition(start_state,bg_state) = 1 - x;
  for(size_t i = first; i <= first + pad_left; i++)
    transition(start_state,i) = x / (pad_left + 1.0);

  // initialize transitions from the background state
  transition(bg_state,start_state) = z;
  transition(bg_state,bg_state) = 1 - y - z;
  for(size_t i = first; i <= first + pad_left; i++)
    transition(bg_state,i) = y / (pad_left + 1.0);

  // initialize transitions from the last of the motif chain states
  for(size_t i = last; i <= last + pad_right; i++) {
    transition(i,start_state) = z;
    transition(i,bg_state) = 1 - y - z;
    for(size_t j = first; j <= first + pad_left; j++)
      transition(i,j) = y / (pad_left + 1.0);
  }
}


void HMM::register_dataset(const Data::Set &data, double class_prior, double motif_p1, double motif_p2)
{
  if(verbosity >= Verbosity::verbose)
    std::cout << "register_data_set(data.path=" << data.path << ( data.is_shuffle ? " shuffle" : " ") << ", sha1=" << data.sha1 << ", class_prior=" << class_prior << ")" << std::endl;
  RegisteredDataSet reg_data({data, class_prior, std::map<size_t, double>()});
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if(is_motif_group(group_idx)) {
      double val = motif_p2;
      if(data.motifs.find(groups[group_idx].name) != end(data.motifs))
        val = motif_p1;
      reg_data.motif_prior[group_idx] = val;
      if(verbosity >= Verbosity::verbose)
        std::cout << "Registering conditional motif prior " << val << " of group " << group_idx << ": " << groups[group_idx].name << " for data set " << data.path << " with sha1 " << data.sha1 << std::endl;
    }
  registered_datasets[data.sha1] = reg_data;
}

