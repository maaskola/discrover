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
 *       Filename:  subhmm.hpp
 *
 *    Description:  Code for a sub model of a hidden Markov model
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef SUBHMM_HPP
#define SUBHMM_HPP

#include <vector>
#include <string>
#include "hmm.hpp"

class SubHMM : public HMM {
public:
  SubHMM(const HMM &hmm, const Training::Range &states);
  size_t n_original_states;
  std::vector<int> reduce;
  std::vector<size_t> lift;
  Training::Range map_down(const Training::Range &range) const;
  Training::Targets map_down(const Training::Targets &targets) const;
  matrix_t lift_emission(const matrix_t &m) const;
  matrix_t lift_transition(const matrix_t &m) const;
  // std::vector<size_t> reduce(const std::vector<size_t> &v) const;
  // std::vector<size_t> lift(const std::vector<size_t> &v) const;
};

#endif
