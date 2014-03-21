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
 *       Filename:  conditional_mutual_information.cpp
 *
 *    Description:  Routines to calculate conditional mutual information
 *
 *        Created:  Wed Mar 21 06:08:55 2014 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include "conditional_mutual_information.hpp"
#include <boost/multi_array.hpp>

using namespace std;

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
double calc_conditional_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity, conditionalMI mode)
{
  const size_t X = 2;
  size_t Y = contrast.sets.size();
  size_t Z = 2;
  if(mode == conditionalMI::motif_and_condition_given_previous_motifs)
    if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info))
      cout << "Conditional mutual information of X=new motif, and Y=contrast, given Z=previous motifs" << endl;
  if(mode == conditionalMI::motif_and_previous_motifs_given_condition) {
    if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info))
      cout << "Conditional mutual information of X=new motif, and Y=previous motifs, given Z=contrast" << endl;
    Y = 2;
    Z = contrast.sets.size();
  }

  typedef boost::multi_array<double, 1> array1_t;
  typedef array1_t::index index;

  typedef boost::multi_array<double, 2> array2_t;
  typedef array2_t::index index;

  typedef boost::multi_array<double, 3> array_t;
  typedef array_t::index index;
  //
  // the joint probability of X, Y, and Z
  array_t p(boost::extents[X][Y][Z]);

  if(mode == conditionalMI::motif_and_condition_given_previous_motifs)
    // fill joint probability table with absolute counts
    for(size_t i = 0; i < contrast.sets.size(); i++) {
      p[0][i][0] = pair_posteriors[i].posterior_both;
      p[0][i][1] = pair_posteriors[i].posterior_first - pair_posteriors[i].posterior_both;
      p[1][i][0] = pair_posteriors[i].posterior_second - pair_posteriors[i].posterior_both;
      p[1][i][1] = pair_posteriors[i].posterior_none;
    }
  else
    // fill joint probability table with absolute counts
    for(size_t i = 0; i < contrast.sets.size(); i++) {
      p[0][0][i] = pair_posteriors[i].posterior_both;
      p[0][1][i] = pair_posteriors[i].posterior_first - pair_posteriors[i].posterior_both;
      p[1][0][i] = pair_posteriors[i].posterior_second - pair_posteriors[i].posterior_both;
      p[1][1][i] = pair_posteriors[i].posterior_none;
    }


  // add pseudo-count
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        p[x][y][z] += ps;

  double marginal = 0;
  // sum over entries to compute marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        marginal += p[x][y][z];

  // normalize joint probability by dividing through marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        p[x][y][z] /= marginal;


  // the joint marginal probability of X and Z
  array2_t pxz(boost::extents[X][Z]);
  // the joint marginal probability of Y and Z
  array2_t pyz(boost::extents[Y][Z]);

  // compute marginal distribution of X and Z by summing over Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pxz[x][z] += p[x][y][z];

  // compute marginal distribution of Y and Z by summing over X
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pyz[y][z] += p[x][y][z];


  // the marginal probability of Z
  array1_t pz(boost::extents[Z]);

  // compute marginal distribution of Z by summing over X and Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        pz[z] += p[x][y][z];

  double mi = 0;
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      for(size_t z = 0; z < Z; z++)
        mi += p[x][y][z] * log( p[x][y][z] / pxz[x][z] / pyz[y][z] * pz[z]);

  mi /= log(2.0);

  mi = max<double>(mi, 0);

  if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info)) {
    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(x,y,z=" << z << "):" << endl;
      for(size_t x = 0; x < X; x++) {
        for(size_t y = 0; y < Y; y++)
          cout << " " << p[x][y][z];
        cout << endl;
      }
    }

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(x,z=" << z << ") =";
      for(size_t x = 0; x < X; x++)
        cout << " " << pxz[x][z];
      cout << endl;
    }

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(y,z=" << z << ") =";
      for(size_t y = 0; y < Y; y++)
        cout << " " << pyz[y][z];
      cout << endl;
    }

    cout << "conditional_mutual_information pz =";
    for(size_t z = 0; z < Z; z++)
      cout << " " << pz[z];
    cout << endl;

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(x,y|z=" << z << "):" << endl;
      for(size_t x = 0; x < X; x++) {
        for(size_t y = 0; y < Y; y++)
          cout << " " << p[x][y][z] / pz[z];
        cout << endl;
      }
    }

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(x|z=" << z << ") =";
      for(size_t x = 0; x < X; x++)
        cout << " " << pxz[x][z] / pz[z];
      cout << endl;
    }

    for(size_t z = 0; z < Z; z++) {
      cout << "conditional_mutual_information p(y|z=" << z << ") =";
      for(size_t y = 0; y < Y; y++)
        cout << " " << pyz[y][z] / pz[z];
      cout << endl;
    }
    cout << "cond-MICO = " << mi << endl;
    cout << endl;
  }

  return(mi);
}

/** The pair mutual information
 * This is the mutual information of the two motifs
 */
double pair_mutual_information(const Data::Contrast &contrast, const HMM::pair_posteriors_t &pair_posteriors, double ps, Verbosity verbosity)
{
  const size_t X = 2;
  const size_t Y = 2;

  typedef boost::multi_array<double, 1> array1_t;
  typedef array1_t::index index;

  typedef boost::multi_array<double, 2> array2_t;
  typedef array2_t::index index;

  // the joint probability of X and Y
  array2_t p(boost::extents[X][Y]);

 // fill joint probability table with absolute counts
  for(size_t i = 0; i < contrast.sets.size(); i++) {
    p[0][0] += pair_posteriors[i].posterior_both;
    p[0][1] += pair_posteriors[i].posterior_first - pair_posteriors[i].posterior_both;
    p[1][0] += pair_posteriors[i].posterior_second - pair_posteriors[i].posterior_both;
    p[1][1] += pair_posteriors[i].posterior_none;
  }

  // add pseudo-count
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      p[x][y] += ps;

  double marginal = 0;
  // sum over entries to compute marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      marginal += p[x][y];

  // normalize joint probability by dividing through marginal
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      p[x][y] /= marginal;


  // the marginal probability of X
  array1_t px(boost::extents[X]);
  // the marginal probability of Y
  array1_t py(boost::extents[Y]);

  // compute marginal distribution of X by summing over Y
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
        px[x] += p[x][y];

  // compute marginal distribution of Y by summing over X
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
        py[y] += p[x][y];

  if(verbosity >= Verbosity::verbose or (verbose_conditional_mico_output and verbosity >= Verbosity::info)) {
    cout << "pair_mutual_information joint =";
    for(size_t x = 0; x < X; x++)
      for(size_t y = 0; y < Y; y++)
        cout << " " << p[x][y];
    cout << endl;
    cout << "pair_mutual_information px =";
    for(size_t x = 0; x < X; x++)
      cout << " " << px[x];
    cout << endl;
    cout << "pair_mutual_information py =";
    for(size_t y = 0; y < Y; y++)
      cout << " " << py[y];
    cout << endl;
  }

  double mi = 0;
  for(size_t x = 0; x < X; x++)
    for(size_t y = 0; y < Y; y++)
      mi += p[x][y] * log( p[x][y] / px[x] / py[y] );

  mi /= log(2.0);

  mi = max<double>(mi, 0);

  return(mi);
}


