/* =====================================================================================
 * Copyright (c) 2014, Jonas Maaskola
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
 *       Filename:  conditional_mutual_information.hpp
 *
 *    Description:  Routines to calculate conditional mutual information
 *
 *        Created:  Wed Mar 21 06:08:55 2014 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef CONDITIONAL_MUTUAL_INFORMATION_HPP
#define CONDITIONAL_MUTUAL_INFORMATION_HPP

#include "hmm.hpp"

enum class conditionalMI {
  motif_and_condition_given_previous_motifs,
  motif_and_previous_motifs_given_condition
};

/** The conditional mutual information
 * This is the reduction in the uncertainty of X due to the knowledge of Y when Z is given:
 * I( X; Y | Z ) = H( X | Z ) - H( X | Y, Z )
 * Symmetrically, it is also the reduction in the uncertainty of Y due to the knowledge of X when Z is given:
 * I( X; Y | Z ) = H( Y | Z ) - H( Y | X, Z )
 *
 * Here X is the first motif, Y is the contrast, and Z is the second (previous) motif(s)
 *
 * See Cover & Thomas 2006 equations (2.60) and (2.61)
 */
double calc_conditional_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity, conditionalMI mode);

/** The pair mutual information
 * This is the mutual information of the two motifs
 */
double pair_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity);

#endif

