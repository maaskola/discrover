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
 *       Filename:  stats.hpp
 *
 *    Description:  Data structure to store counts
 *
 *        Created:  Thu Jun 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef  STATS_HPP
#define  STATS_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef double fp_t;
typedef boost::numeric::ublas::matrix<fp_t> matrix_t;
typedef boost::numeric::ublas::vector<fp_t> vector_t;
typedef boost::numeric::ublas::matrix<size_t> count_matrix_t;
typedef boost::numeric::ublas::vector<size_t> count_vector_t;

typedef boost::numeric::ublas::zero_matrix<fp_t> zero_matrix_t;
typedef boost::numeric::ublas::zero_matrix<size_t> zero_count_matrix_t;
typedef boost::numeric::ublas::zero_vector<fp_t> zero_vector_t;
typedef boost::numeric::ublas::zero_vector<size_t> zero_count_vector_t;

typedef boost::numeric::ublas::scalar_matrix<fp_t> scalar_matrix_t;
typedef boost::numeric::ublas::scalar_matrix<size_t> scalar_count_matrix_t;
typedef boost::numeric::ublas::scalar_vector<fp_t> scalar_vector_t;
typedef boost::numeric::ublas::scalar_vector<size_t> scalar_count_vector_t;

namespace Seeding {
  namespace Stats {
    // typedef std::vector<size_t> Stats;
    // typedef std::vector<size_t> OccurrenceCounts;
    typedef count_vector_t OccurrenceCounts;
    typedef matrix_t OccurrenceTable;
    // struct Contrast{
    //   Stats signal, control;
    // };
  }
}

#endif   /* ----- #ifndef STATS_HPP ----- */

