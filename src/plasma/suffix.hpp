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
 *       Filename:  suffix.hpp
 *
 *    Description:  Generic code for Suffix array based alignment
 *
 *        Created:  Tue Aug 08 22:23:44 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef SUFFIX_HPP
#define SUFFIX_HPP

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <stack>
#include <numeric>
#include <algorithm>
#include "../timer.hpp"
#include "../aux.hpp"
#include "construction.hpp"
#include "../verbosity.hpp"

/** A direct but slow (n log(n)) algorithm to construct the suffix array */
template <class idx_t, class Iter>
std::vector<idx_t> gen_suffix_array_slow(Iter begin, Iter end,
                                         Verbosity verbosity) {
  Timer timer;
  std::vector<idx_t> sa(std::distance(begin, end));
  std::iota(sa.begin(), sa.end(), 0);
  auto cmp = [begin, end](idx_t a, idx_t b) {
    return lexicographical_compare(begin + a, end, begin + b, end);
  };
  std::sort(sa.begin(), sa.end(), cmp);
  double time = timer.tock();
  if (verbosity >= Verbosity::verbose)
    std::cerr << "Building SA took " + to_pretty_string(time) + " µs." << std::endl;
  return sa;
}

template <typename X>
X plusOne(X i) {
  return ++i;
}

/** A faster (linear-time) algorithm to construct the suffix array */
template <bool shift, class idx_t, class Iter>
std::vector<idx_t> gen_suffix_array(Iter begin, const Iter end,
                                    Verbosity verbosity) {
  idx_t K = 0;
  if (begin != end)
    K = *std::max_element(begin, end);
  const idx_t n = std::distance(begin, end);
  Timer timer;

  // we need to copy the data and make sure that it is padded with three zeros
  // for the DC3 algorithm to work
  // also, zeros are not allowed, so optionally we increment the data
  std::vector<typename Iter::value_type> v(n + 3, 0);
  if (shift)
    std::transform(begin, end, v.begin(), plusOne<typename Iter::value_type>);
  else
    std::copy(begin, end, v.begin());

  std::vector<idx_t> sa(n + 3, 0);
  suffixArray(v.begin(), v.end(), sa, n, K + (shift ? 1 : 0));
  double time = timer.tock();
  if (verbosity >= Verbosity::verbose)
    std::cerr << "Building SA took " + to_pretty_string(time) + " µs." << std::endl;
  return sa;
}

/** A simple but slow algorithm to construct the LCP array */
template <class lcp_t, class idx_t, class Iter>
std::vector<lcp_t> gen_lcp_slow(Iter begin, Iter end,
                                const std::vector<idx_t> &sa,
                                Verbosity verbosity) {
  Timer timer;
  idx_t n = sa.size();
  std::vector<lcp_t> lcp(n, 0);
  for (idx_t i = 1; i < n; i++) {
    Iter a = begin + sa[i - 1];
    Iter b = begin + sa[i];
    while (a != end and b != end and *(a++) == *(b++))
      lcp[i]++;
  }
  double time = timer.tock();
  if (verbosity >= Verbosity::verbose)
    std::cerr << "Building LCP took " + to_pretty_string (time) + " µs." << std::endl;
  return lcp;
}

/*  Pseudocode of linear time LCP array creation due to Kasai et al
 *  T. Kasai, G. Lee, H. Arimura, S. Arikawa, and K. Park.
 *  Linear-time longest-common-prefix computation in suffix arrays and its applications.
 *  In Proc. CPM, volume 2089 of LNCS, pages 181–192. Springer, 2001.
 *
 * Algorithm GetHeight
 * input: a text A and its suffix array Pos
 * output: the LCP table Height
 *
 * for i:=1 to n:
 *   Rank[Pos[i]] := i
 *  h:=0
 *  for i:=1 to n:
 *    if Rank[i] > 1:
 *      k := Pos[Rank[i]-1]         // note: in the manuscript it says "k" here; this is likely a mistake and it needs to be "j" instead!
 *      while A[i+h] = A[j+h]:
 *        h := h+1
 *      Height[Rank[i]] := h
 *      if h > 0:
 *        h := h-1
 */
template <class lcp_t, class idx_t, class Iter>
std::vector<lcp_t> gen_lcp(Iter begin, Iter end, const std::vector<idx_t> &sa,
                           Verbosity verbosity) {
  Timer timer;
  idx_t n = std::distance(begin, end);
  std::vector<idx_t> lcp(n), rank(n);
  for (idx_t i = 0; i < n; i++)
    rank[sa[i]] = i;
  idx_t h = 0;
  for (idx_t i = 0; i < n; i++)
    if (rank[i] > 0) {
      idx_t j = sa[rank[i] - 1];
      while (*(begin + i + h) == *(begin + j + h))
        h++;
      lcp[rank[i]] = h;
      if (h > 0)
        h--;
    }
  double time = timer.tock();
  if (verbosity >= Verbosity::verbose)
    std::cerr << "Building LCP took " + to_pretty_string(time) + " µs." << std::endl;
  return lcp;
}

// TODO use more efficient algorithm (although so far this is not a bottle-neck)
// TODO another thought should be spent on the definition of the jmp pointer;
// right now it's defined to be the next suffix with lcp <= to the current;
// perhaps it might be helpful to use strictly less instead.
template <class idx_t, class lcp_t>
std::vector<idx_t> gen_jmp(const std::vector<lcp_t> &lcp, Verbosity verbosity) {
  Timer timer;
  idx_t n = lcp.size();
  std::vector<idx_t> jmp(n);
  for (idx_t i = 0; i < n; i++) {
    jmp[i] = i + 1;
    size_t j = i + 1;
    while (j < n and lcp[j++] > lcp[i])
      jmp[i]++;
  }
  double time = timer.tock();
  if (verbosity >= Verbosity::verbose)
    std::cerr << "Building JMP took " + to_pretty_string(time) + " µs." << std::endl;
  return jmp;
}

template <class idx_t, class lcp_t>
std::vector<idx_t> gen_jmp_var(const std::vector<lcp_t> &lcp) {
  idx_t n = lcp.size();
  std::vector<idx_t> jmp(n);

  std::stack<std::pair<idx_t, lcp_t>> s;
  s.push(std::make_pair(0, 0));
  for (idx_t i = 1; i < n; i++) {
    while (lcp[i] < s.top().second) {
      jmp[s.top().first] = i;
      s.pop();
    }
    s.push(std::make_pair(i, lcp[i]));
  }
  while (not s.empty()) {
    jmp[s.top().first] = n;
    s.pop();
  }
  return jmp;
}

template <class lcp_t, class idx_t, class Iter,
          typename Cmp = std::equal_to<typename Iter::value_type>>
std::vector<idx_t> match(Iter qbegin, Iter qend, Iter begin, Iter end,
                         const std::vector<idx_t> &sa,
                         const std::vector<lcp_t> &lcp,
                         const std::vector<idx_t> &jmp, Cmp cmp = Cmp()) {
  Timer timer;
  const bool do_debug = false;
  std::vector<idx_t> hits;

  const size_t qsize = std::distance(qbegin, qend);
  std::stack<idx_t> jmp_stack;
  const idx_t n = std::distance(begin, end);

  idx_t i = 0;
  while (i < n) {
    if (do_debug)
      std::cout << "i = " << i << std::endl;
    // invariant: lcp[i] > length of current match
    // consequence: there is always at least one character to match
    jmp_stack.push(jmp[i]);
    lcp_t j = lcp[i];
    auto q = qbegin + j;
    auto t = begin + sa[i] + j;
    if (do_debug)
      std::cout << "cmp = " << *t << " " << *q << std::endl;
    while (q != qend and t != end and cmp(*t, *q)) {
      if (do_debug)
        std::cout << "matched = " << *t << " " << *q << std::endl;
      q++;
      t++;
      if (do_debug and q != qend and t != end)
        std::cout << "cmp = " << *t << " " << *q << std::endl;
    }
    if (q == qend) {  // found a match
      if (do_debug)
        std::cout << "adding hit: " << sa[i] << std::endl;
      hits.push_back(sa[i]);
      i++;
      while (i < n and lcp[i] >= qsize) {
        idx_t j = i;
        i = jmp[i];
        for (; j < i; j++) {
          hits.push_back(sa[j]);
          if (do_debug)
            std::cout << "adding subsequent hit: " << sa[j] << std::endl;
        }
      }
    } else if (q == qbegin + j) {  // mismatch on first comparison
      i = jmp_stack.top();
      jmp_stack.pop();
      if (do_debug)
        std::cout << "Jump!" << std::endl;
    } else {  // mismatch on later position
      i++;
      size_t k = std::distance(qbegin, q);
      if (do_debug)
        std::cout << "i = " << i << " k = " << k << " lcp = " << lcp[i]
                  << std::endl;
      while (i < n and lcp[i] > k) {
        i = jmp[i];
        if (do_debug)
          std::cout << "Bounce!" << std::endl;
        if (do_debug)
          std::cout << "i = " << i << " k = " << k << " lcp = " << lcp[i]
                    << std::endl;
      }
    }
  }
  double time = timer.tock();
  if (do_debug) {
    std::cerr << "Getting matches took " + to_pretty_string(time) + " µs." << std::endl;
    std::cerr << "Got " << hits.size() << " results." << std::endl;
  }
  return hits;
}

#endif
