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
 *       Filename:  hmm.cpp
 *
 *    Description:  Constructors of the HMM class 
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <fstream>
#include <boost/filesystem.hpp>
#include "aux.hpp"
#include "hmm.hpp"

const size_t HMM::start_state;
const size_t HMM::bg_state;

HMM::HMM(const std::string &path, Verbosity verbosity_, double pseudo_count) :
  n_emissions(alphabet_size),
  verbosity(verbosity_),
  store_intermediate(false),
  last_state(0),
  n_states(0),
  pseudo_count(pseudo_count),
  groups(),
  group_ids(),
  order(),
  order_offset(),
  max_order(0),
  transition(),
  emission(),
  pred(),
  succ(),
  all_range(),
  constitutive_range(),
  emitting_range(),
  registered_datasets()
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Called HMM constructor 1." << std::endl;
  if(not boost::filesystem::exists(path)) {
    std::cout << "Error: HMM parameter file " << path << " does not exist. Aborting." << std::endl;
    exit(-1);
  }
  std::ifstream ifs(path.c_str());
  if(not ifs.good()) {
    std::cout << "Error: can't read from parameter file " << path << ". Aborting." << std::endl;
    exit(-1);
  }
  try {
    deserialize(ifs);
  } catch (...) {
    std::cout << "Error during loading parameters from parameter file " << path << ". Aborting." << std::endl;
    exit(-1);
  }
  initialize_order_offsets();
  finalize_initialization();
};

HMM::HMM(const HMM &hmm, bool copy_deep) :
  n_emissions(hmm.n_emissions),
  verbosity(hmm.verbosity),
  store_intermediate(hmm.store_intermediate),
  last_state(hmm.last_state),
  n_states(hmm.n_states),
  pseudo_count(hmm.pseudo_count),
  groups(hmm.groups),
  group_ids(hmm.group_ids),
  order(hmm.order),
  order_offset(hmm.order_offset),
  max_order(hmm.max_order),
  transition(hmm.transition),
  emission(hmm.emission),
  pred(hmm.pred),
  succ(hmm.succ),
  all_range(hmm.all_range),
  constitutive_range(hmm.constitutive_range),
  emitting_range(hmm.emitting_range),
  registered_datasets(hmm.registered_datasets)
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Called HMM constructor 2." << std::endl;
  initialize_order_offsets();
  finalize_initialization();
};

size_t calc_num_emissions(size_t asize, size_t order)
{
  size_t x = 0;
  for(size_t i = 1; i <= order+1; i++)
    x += ipow(asize, i);
  return(x);
}

HMM::HMM(size_t bg_order, Verbosity verbosity_, double pseudo_count_) :
  n_emissions(calc_num_emissions(alphabet_size, bg_order)),
  verbosity(verbosity_),
  store_intermediate(false),
  last_state(1),
  n_states(2), // for the start state and background
  pseudo_count(pseudo_count_),
  groups(),
  group_ids(),
  order(std::vector<size_t>(n_states, 0)),
  order_offset(std::vector<size_t>(n_states)),
  max_order(bg_order),
  transition(zero_matrix(n_states, n_states)),
  emission(zero_matrix(n_states, n_emissions)),
  pred(),
  succ(),
  all_range(),
  constitutive_range(),
  emitting_range(),
  registered_datasets()
{
  if(verbosity >= Verbosity::debug)
    std::cout << "Called HMM constructor 3." << std::endl
      << "bg_order = " << bg_order << std::endl;
  Group special({Group::Kind::Special, "Special", {start_state}});
  groups.push_back(special);
  group_ids.push_back(0);

  Group background({Group::Kind::Background, "Background", {bg_state}});
  groups.push_back(background);
  group_ids.push_back(1);
  order[bg_state] = bg_order; // for the start and end state
  initialize_order_offsets();
  initialize_transitions();
  initialize_bg_transitions();

  initialize_emissions(bg_order);

  normalize_emission(emission);
  normalize_transition(transition);
  finalize_initialization();
  if(verbosity >= Verbosity::debug)
    std::cout << "Constructed HMM: " << *this << std::endl;
};

size_t HMM::add_motif(const matrix_t &e, double exp_seq_len, double lambda, const std::string &name, const std::vector<size_t> &insertions)
{
  if(verbosity >= Verbosity::verbose)
    std::cout << "Adding motif " << name << ": " << e << std::endl;

  size_t n_insertions = insertions.size();
  size_t n = e.size1() + n_insertions;

  HMM previous(*this);

  size_t motif_idx = groups.size();

  size_t first = n_states;
  size_t last = first + n - 1;
  Group motif = {Group::Kind::Motif, name, {}};
  for(size_t i = first ; i <= last; i++) {
    motif.states.push_back(i);
    group_ids.push_back(motif_idx);
  }
  groups.push_back(motif);

  last_state += n;
  n_states += n;

  transition = zero_matrix(n_states, n_states);
  emission = zero_matrix(n_states, n_emissions);

  for(size_t i = 0; i < previous.transition.size1(); i++)
    for(size_t j = 0; j < n_emissions; j++)
      emission(i,j) = previous.emission(i,j);

  for(size_t i = 0; i < previous.transition.size1(); i++)
    for(size_t j = 0; j < previous.transition.size2(); j++)
      transition(i,j) = previous.transition(i,j);

  for(size_t i = 0; i < n; i++)
    order.push_back(0); // TODO: make configurable

  set_motif_emissions(e, first, n_insertions);

  initialize_transitions_to_and_from_chain(n, exp_seq_len, lambda, first, last);

  // initialize transitions between successive chain states
  for(size_t i = first; i < last; i++)
    transition(i,i+1) = 1.0;

  normalize_transition(transition);
  normalize_emission(emission);

  if(verbosity >= Verbosity::verbose) {
    std::cout << "Insertions:";
    for(auto &x: insertions)
      std::cout << " " << x;
    std::cout << std::endl;
  }
  if(n_insertions > 0) {
    for(auto &i: insertions)
      if(i == 0 or i > e.size1()) {
        // insertions are currently only allowed between motif states
        // TODO: generalize
        std::cout << "Error: insertion positions must be within the motif.";
        exit(-1);
      }


    const double insert_transition_probability = 0.5;

    size_t first_insertion = last - n_insertions;
    size_t second_to_last_regular = first_insertion - 1;

    if(verbosity >= Verbosity::debug)
      std::cout << transition << std::endl;

    // set to zero the transition probabilities of insertion states
    for(size_t i = first_insertion; i < last; i++)
      for(size_t j = 0; j < n_states; j++)
        transition(i,j) = 0;

    if(verbosity >= Verbosity::debug)
      std::cout << transition << std::endl;

    // fix transition of second to last chain state
    for(size_t j = 0; j < n_states; j++)
      transition(second_to_last_regular, j) = 0;
    transition(second_to_last_regular, last) = 1;

    if(verbosity >= Verbosity::debug)
      std::cout << transition << std::endl;

    // fix transitions of insertions
    for(size_t i_ = 0; i_ < n_insertions; i_++) {
      size_t i = n_insertions - 1 - i_;
      for(size_t j = 0; j < n_states; j++) {
        transition(first_insertion+i, j) = transition(insertions[i]+1,j);
        transition(insertions[i]+1,j) *= (1 - insert_transition_probability);
      }
      transition(insertions[i]+1,first_insertion+i) = insert_transition_probability;
    }
    if(verbosity >= Verbosity::debug)
      std::cout << transition << std::endl;
  }

  normalize_transition(transition);
  normalize_emission(emission);

  finalize_initialization();

  return(motif_idx);
}

void set_emissions(size_t i, const std::vector<size_t> &these, matrix_t &m, double x)
{
  double alpha = (1.0 - x * (4-these.size())) / these.size();
  for(size_t j = 0; j < 4; j++)
    m(i,j) = x;
  for(size_t j = 4; j < m.size2(); j++)
    m(i,j) = 0;
  for(auto j : these)
    m(i,j) = alpha;
  // normalize to be safe
  double z = 0;
  for(size_t j = 0; j < m.size2(); j++)
    z += m(i,j);
  for(size_t j = 0; j < m.size2(); j++)
    m(i,j) /= z;
}

size_t HMM::add_motif(const std::string &seq, double alpha, double exp_seq_len, double lambda, const std::string &name, const std::vector<size_t> &insertions)
{
  if(verbosity >= Verbosity::verbose)
    std::cout << "Adding motif for " << seq << "." << std::endl;
  matrix_t e = zero_matrix(seq.size(), n_emissions);
  for(size_t i = 0; i < seq.size(); i++) {
    char c = tolower(seq[i]);
    std::vector<size_t> these;
    switch(c) {
      case 'a':
        these = {0};
        break;
      case 'c':
        these = {1};
        break;
      case 'g':
        these = {2};
        break;
      case 't':
      case 'u':
        these = {3};
        break;
      case 'r':
        these = {0,2};
        break;
      case 'y':
        these = {1,3};
        break;
      case 's':
        these = {1,2};
        break;
      case 'w':
        these = {0,3};
        break;
      case 'k':
        these = {2,3};
        break;
      case 'm':
        these = {0,1};
        break;
      case 'b':
        these = {1,2,3};
        break;
      case 'd':
        these = {0,2,3};
        break;
      case 'h':
        these = {0,1,3};
        break;
      case 'v':
        these = {0,1,2};
        break;
      case 'n':
        these = {0,1,2,3};
        break;
      default:
        std::cout << "Error: unknown encoding!" << std::endl;
        throw("Unknown encoding");
    }
    set_emissions(i, these, e, alpha);
  }
  size_t motif_idx = add_motif(e, exp_seq_len, lambda, name, insertions);
  return(motif_idx);
}

void HMM::add_motifs(const HMM &hmm) {
  for(auto &group: hmm.groups)
    if(group.kind == Group::Kind::Motif) {
      size_t n_motif_states = group.states.size();
      matrix_t new_transition = zero_matrix(n_states + n_motif_states, n_states + n_motif_states);
      matrix_t new_emission = zero_matrix(n_states + n_motif_states, n_emissions);
      for(size_t i = 0; i < n_states; i++)
        for(size_t j = 0; j < n_states; j++)
          new_transition(i,j) = transition(i,j);
      for(size_t i = n_states; i < n_states + n_motif_states - 1; i++)
        new_transition(i, i+1) = 1;
      new_transition(bg_state, n_states) = hmm.transition(hmm.bg_state, *group.states.begin());
      new_transition(bg_state, bg_state) -= new_transition(bg_state, n_states);
      new_transition(n_states + n_motif_states - 1, bg_state) = hmm.transition(*group.states.rbegin(), hmm.bg_state);
      new_transition(n_states + n_motif_states - 1, n_states) = 1 - new_transition(n_states + n_motif_states - 1, bg_state);

      for(size_t i = 0; i < n_states; i++)
        for(size_t j = 0; j < n_emissions; j++)
          new_emission(i,j) = emission(i,j);
      for(size_t i = 0; i < n_motif_states; i++)
        for(size_t j = 0; j < n_emissions; j++)
          new_emission(n_states + i, j) = hmm.emission(group.states[i], j);

      size_t group_idx = groups.size();
      for(size_t i = 0; i < n_motif_states; i++) {
        group_ids.push_back(group_idx);
        order.push_back(hmm.order[group.states[i]]);
      }

      Group new_group = group;
      for(size_t i = 0; i < n_motif_states; i++)
        new_group.states[i] = n_states + i;
      groups.push_back(new_group);

      n_states += n_motif_states;
      last_state += n_motif_states;

      emission = new_emission;
      transition = new_transition;
    }
  finalize_initialization();
}

