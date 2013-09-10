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
 *       Filename:  aux.hpp
 *
 *    Description:  Auxiliary routines
 *
 *        Created:  Tue Aug 8 22:23:44 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef AUX_HPP
#define AUX_HPP

#include <string>
#include <vector>

template <typename Iter> void range_tolower(Iter beg, Iter end) {
  for(Iter iter = beg; iter != end; ++iter) {
    *iter = std::tolower(*iter);
  }
}

std::string string_tolower(const std::string & str);

template <class T> T hibit_orig(T n) {
  n |= (n >>   1);
  n |= (n >>   2);
  n |= (n >>   4);
  n |= (n >>   8);
  n |= (n >>  16);
  n |= (n >>  32);
  return n - (n >> 1);
};

template <class T> T hibit(T n) {
  if(n == 0)
    return(0);
  size_t bits = 1;
  while((n = (n >> 4)) > 0)
    bits++;
  return bits;
};

/** Parse a comma separated list of ranges.
 * Example: "1,2,5-7,2,5,30"  will yield {1, 2, 5, 6, 7, 2, 5, 30}. */
std::vector<std::string> tokenize(const std::string &s, const std::string &delim);
std::vector<size_t> parse_list(const std::string &s);

std::string sha1hash(const std::string &s);

#endif   /* ----- #ifndef AUX_HPP  ----- */

