/*
 * =====================================================================================
 *
 *       Filename:  stats.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */


#ifndef  STATS_HPP
#define  STATS_HPP

//#include <vector>
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

namespace Plasma {
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

