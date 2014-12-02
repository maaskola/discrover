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
 *       Filename:  conditional.hpp
 *
 *    Description:  Find the best motif occurrences
 *
 *        Created:  Tue Dec 2 23:05:23 2014 +0100
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef CONDITIONAL_DECODER_HPP
#define CONDITIONAL_DECODER_HPP

#include <iostream>
#include "hmm.hpp"

class ConditionalDecoder {
  HMM hmm;
  std::vector<std::pair<std::string, std::list<matrix_t>>> emission_matrices;

  public:
  ConditionalDecoder(const HMM &hmm_);
  void decode(std::ostream &os, const Data::Seq &seq) const;

  protected:
};

#endif
