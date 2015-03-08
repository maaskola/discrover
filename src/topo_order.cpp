/* =====================================================================================
 * Copyright (c) 2015, Jonas Maaskola
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
 *       Filename:  topo_order.cpp
 *
 *    Description:  Compute topological ordering
 *
 *        Created:  Sun Mar 8 16:11:23 2015 +0100
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "topo_order.hpp"
#include <map>
#include <set>

using namespace std;

// Find those states reachable from the set of states given as argument.
// NOTE: There's a slight hack in that transitions from the set of final
// states to the set of initial states are not considered.
map<size_t, set<size_t>> reachable_states(const matrix_t &transition,
                                          const vector<size_t> &states,
                                          const set<size_t> &initial,
                                          const set<size_t> &final) {
  map<size_t, set<size_t>> r;
  for (auto &i : states)
    for (auto &j : states)
      if (transition(i, j) > 0)
        if (final.find(i) == end(final) or initial.find(j) == end(initial))
          r[i].insert(j);
  size_t n = states.size();
  while (n--) {
    for (auto &p : r) {
      set<size_t> y;
      for (auto &z : p.second) {
        y.insert(z);
        for (auto &a : r[z])
          y.insert(a);
      }
      p.second = y;
    }
  }
  return r;
}

/** Find states reachable from a given state. */
set<size_t> reachable_from_state(const matrix_t &transition,
                                 const size_t &state) {
  set<size_t> r;
  for (size_t i = 0; i < transition.size2(); i++)
    if (transition(state, i) > 0)
      r.insert(i);
  return r;
}

/** Determine initial motif states.
 *  These are those with connections from to the background state.
 */
set<size_t> initial_states(const matrix_t &transition,
                           const vector<size_t> &states) {
  const size_t bg_state = 1;  // HMM::bg_state is protected
  auto reachable_from_bg = reachable_from_state(transition, bg_state);
  set<size_t> initial;
  for (auto state : states)
    if (reachable_from_bg.find(state) != end(reachable_from_bg))
      initial.insert(state);
  return initial;
}

/** Determine final motif states.
 *  These are those with connections back to the background state.
 */
set<size_t> final_states(const matrix_t &transition,
                         const vector<size_t> &states) {
  const size_t bg_state = 1;  // HMM::bg_state is protected
  set<size_t> final;
  for (auto state : states)
    if (transition(state, bg_state) > 0)
      final.insert(state);
  return final;
}

vector<size_t> topological_order(const matrix_t &transition,
                                 const vector<size_t> &states) {
  auto initial = initial_states(transition, states);
  auto final = final_states(transition, states);
  auto reachable = reachable_states(transition, states, initial, final);
  vector<size_t> order = states;
  sort(begin(order), end(order), [&reachable](size_t x, size_t y) {
    if (x == y)
      return true;
    else if (reachable[x].find(y) != end(reachable[x]))
      return true;
    else
      return false;
  });

  return order;
}
