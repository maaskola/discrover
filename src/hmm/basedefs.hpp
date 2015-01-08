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
 *       Filename:  basedefs.hpp
 *
 *    Description:  Basic classes and type definitions for sets of sequence sets
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef BASEDEFS_HPP
#define BASEDEFS_HPP

#include <string>
#include "../plasma/data.hpp"
#include "trainingmode.hpp"
#include "sequence.hpp"

namespace Data {
using Seq = Fasta::IEntry;
using Set = Basic::Set<Seq>;
using Contrast = Basic::Contrast<Set>;
using Collection = Basic::Collection<Contrast>;

using Seqs = std::vector<Seq>;
}

void prepare_cross_validation(const Data::Collection &col,
                              Data::Collection &training_data,
                              Data::Collection &test_data,
                              double cross_validation_freq,
                              Verbosity verbosity);
void prepare_cross_validation(const Data::Contrast &contrast,
                              Data::Contrast &training_data,
                              Data::Contrast &test_data,
                              double cross_validation_freq,
                              Verbosity verbosity);

#endif
