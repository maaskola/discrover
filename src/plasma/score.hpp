/*
 * =====================================================================================
 *
 *       Filename:  score.hpp
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


#ifndef  SCORE_HPP
#define  SCORE_HPP

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "plasma_stats.hpp"
#include "motif.hpp"
#include "options.hpp"
#include "results.hpp"
// #include "../GitSHA1.hpp"

typedef double fp_t;
typedef boost::numeric::ublas::matrix<fp_t> matrix_t;
typedef boost::numeric::ublas::vector<fp_t> vector_t;


double compute_mutual_information_variance(const Plasma::Stats::OccurrenceTable &m_, double pseudo_count, bool normalize);
double compute_mutual_information(const Plasma::Stats::OccurrenceTable &counts, double pseudo_count=0, bool normalize=false, bool do_correction=false);
double compute_mutual_information(double a, double b, double c, double d);

double compute_mcc(double a, double b, double c, double d);

double compute_delta_frequency(double a, double b, double c, double d);

/** Compute the score
 * If measure is Measure::Undefined then the measure is used that is a member of options_t options
 **/
double compute_score(
    const Plasma::DataCollection &collection,
    const Plasma::Result &result,
    const Plasma::options_t &options,
    Measures::Discrete::Measure measure=Measures::Discrete::Measure::Undefined,
    bool do_correction=false
  );
double compute_score(
    const Plasma::DataCollection &collection,
    const Plasma::Stats::OccurrenceCounts &counts,
    const Plasma::options_t &options,
    const Plasma::Objective &objective,
    size_t length,
    size_t degeneracy,
    Measures::Discrete::Measure measure=Measures::Discrete::Measure::Undefined,
    bool do_correction=false);
double compute_score(
    const Plasma::DataSeries &data_series,
    const Plasma::Stats::OccurrenceCounts &contrast,
    const Plasma::options_t &options,
    Measures::Discrete::Measure measure,
    size_t length,
    size_t degeneracy,
    const std::string &motif_name="",
    bool do_correction=false);

double approximate_score(const std::string &motif, const Plasma::hash_map_t &counts, const Plasma::options_t &options);

#endif   /* ----- #ifndef SCORE_HPP  ----- */

