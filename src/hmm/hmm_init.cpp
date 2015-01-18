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
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "hmm.hpp"

using namespace std;

void HMM::finalize_initialization() {
  initialize_pred_succ();
  check_consistency();
}

void HMM::initialize_pred_succ() {
  pred = vector<vector<size_t>>();
  succ = vector<vector<size_t>>();
  for (size_t i = 0; i < n_states; i++) {
    pred.push_back(vector<size_t>());
    succ.push_back(vector<size_t>());
  }
  for (size_t i = 0; i < n_states; i++)
    for (size_t j = 0; j < n_states; j++)
      if (transition(i, j) > 0) {
        pred[j].push_back(i);
        succ[i].push_back(j);
      }
};

/** Initialize the emission matrix */
void HMM::initialize_emissions() {
  // Start state does not emit
  for (size_t i = 0; i < n_emissions; i++)
    emission(start_state, i) = 0;

  // All other states are uniform over all emissions
  for (size_t i = bg_state; i < n_states; i++)
    for (size_t j = 0; j < n_emissions; j++)
      emission(i, j) = 1.0 / n_emissions;
};

/** Initialize the transition matrix to zero */
void HMM::initialize_transitions() {
  transition = zero_matrix(n_states, n_states);
}

void HMM::initialize_bg_transitions() {
  transition(start_state, bg_state) = 1.0;
  transition(bg_state, start_state) = 1.0;
  transition(bg_state, bg_state) = 1.0;
};

/** Initialize the transition probabilities to and from the motif chain such
  * that lambda of the sequences have a motif occurrence, where the expected
  * sequence length is l.
  */
void HMM::initialize_transitions_to_and_from_chain(size_t w, double l,
                                                   double lambda, size_t first,
                                                   size_t last, size_t pad_left,
                                                   size_t pad_right) {
  const double x = lambda / (l - w + 1);
  const double y = lambda / (l - w + 1) * (l - w) / (l - w * lambda);
  const double z = (1 - lambda / (l - w + 1)) / (l - w * lambda);

  if (verbosity >= Verbosity::debug)
    cout << "l = " << l << endl << "w = " << w << endl << "x = " << x << endl
         << "y = " << y << endl << "z = " << z << endl;

  // initialize transitions from the start state
  transition(start_state, start_state)
      = 0;  // sequences are at least one position long
  transition(start_state, bg_state) = 1 - x;
  for (size_t i = first; i <= first + pad_left; i++)
    transition(start_state, i) = x / (pad_left + 1.0);

  // initialize transitions from the background state
  transition(bg_state, start_state) = z;
  transition(bg_state, bg_state) = 1 - y - z;
  for (size_t i = first; i <= first + pad_left; i++)
    transition(bg_state, i) = y / (pad_left + 1.0);

  // initialize transitions from the last of the motif chain states
  for (size_t i = last; i <= last + pad_right; i++) {
    transition(i, start_state) = z;
    transition(i, bg_state) = 1 - y - z;
    for (size_t j = first; j <= first + pad_left; j++)
      transition(i, j) = y / (pad_left + 1.0);
  }
}
