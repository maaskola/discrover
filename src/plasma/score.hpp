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
 *       Filename:  score.hpp
 *
 *    Description:  Discriminative objective functions
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef SCORE_HPP
#define SCORE_HPP

#include "plasma_stats.hpp"
#include "options.hpp"
#include "results.hpp"
#include "data.hpp"

double compute_mutual_information_variance(const matrix_t &m_, double pseudo_count, bool normalize);
double compute_mutual_information(const matrix_t &counts, double pseudo_count=0, bool normalize=false, bool do_correction=false);
double compute_mutual_information(double a, double b, double c, double d);

double compute_mcc(double a, double b, double c, double d);

double compute_delta_frequency(double a, double b, double c, double d);

/** Compute the score
 * If measure is Measure::Undefined then the measure is used that is a member of Options options
 **/
double compute_score(
    const Seeding::Collection &collection,
    const Seeding::Result &result,
    const Seeding::Options &options,
    Measures::Discrete::Measure measure=Measures::Discrete::Measure::Undefined,
    bool do_correction=false
  );
double compute_score(
    const Seeding::Collection &collection,
    const count_vector_t &counts,
    const Seeding::Options &options,
    const Seeding::Objective &objective,
    size_t length,
    size_t degeneracy,
    Measures::Discrete::Measure measure=Measures::Discrete::Measure::Undefined,
    bool do_correction=false);
double compute_score(
    const Seeding::Contrast &contrast,
    const count_vector_t &counts,
    const Seeding::Options &options,
    Measures::Discrete::Measure measure,
    size_t length,
    size_t degeneracy,
    const std::string &motif_name="",
    bool do_correction=false);

double approximate_score(const std::string &motif, const Seeding::hash_map_t &counts, const Seeding::Options &options);

#endif   /* ----- #ifndef SCORE_HPP  ----- */

