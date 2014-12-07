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
 *       Filename:  code.hpp
 *
 *    Description:  IUPAC nucleic acid codes
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef  CODE_HPP
#define  CODE_HPP

#include <string>
#include <vector>

namespace Seeding {
  const std::string Symbol = "-acmgrsvtwyhkdbn";
  std::string iupac2regex(const std::string &s);
  bool iupac_included(char r, char q);
};

using symbol_t = uint8_t;
using seq_type = std::vector<symbol_t>;

/** Encode nucleic acid sequences
 * Takes strings and produces seq_type
 * IUPAC characters (whether upper or lower case) are faithfully encoded
 * All other characters are represented by non-matching symbols
 **/
seq_type encode(const std::string &seq);
/** Decode nucleic acid sequences, producing lower case strings
 **/
std::string decode(const seq_type &seq);
/** Extends seq_type nucleic acid sequences by sequences given in strings
 * Optionally, IUPAC wildcards are faithfully encoded.
 * Otherwise, they and all non-IUPAC symbols are encoded by non-matching symbols
 */
void add_sequence(seq_type &s, const std::string &seq, bool allow_iupac_wildcards);

bool pure_nucleotide(symbol_t s);
bool degenerate_nucleotide(symbol_t s);

seq_type iupac_reverse_complement(const seq_type &s);

#endif /* ----- #ifndef CODE_HPP ----- */
