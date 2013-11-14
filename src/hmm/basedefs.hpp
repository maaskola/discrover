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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
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
  typedef Fasta::IEntry Seq;
  typedef Basic::Set<Seq> Set;
  typedef Basic::Series<Set> Series;
  typedef Basic::Collection<Series> Collection;

  typedef std::vector<Seq> Seqs;
}

namespace Training {
  struct State {
    int center;
    std::vector<std::vector<double>> scores;
    State(size_t n);
  };
}

void prepare_cross_validation(const Data::Collection &data_sets, Data::Collection &training_data, Data::Collection &test_data, double cross_validation_freq, Verbosity verbosity);
void prepare_cross_validation(const Data::Series &data, Data::Series &training_data, Data::Series &test_data, double cross_validation_freq, Verbosity verbosity);

#endif

