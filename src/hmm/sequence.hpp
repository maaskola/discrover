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
 *       Filename:  sequence.hpp
 *
 *    Description:  Typedefs and routines for nucleic acid sequences
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <string>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include "../verbosity.hpp"

typedef unsigned char alphabet_idx_t;
typedef boost::numeric::ublas::vector<alphabet_idx_t> seq_t;

const size_t empty_symbol = 5;
seq_t string2seq(const std::string &s, int n_enc=-1);
std::string seq2string(const seq_t &s);
seq_t random_seq(size_t n, size_t alphabet_size);

std::vector<std::string> extract_seq_ids(const std::string &path, size_t nseq, Verbosity verbosity);


#endif

