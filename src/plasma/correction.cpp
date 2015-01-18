/*
 * =====================================================================================
 *
 *       Filename:  score.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "score.hpp"

using namespace std;

matrix_t compute_zeta_table() {
  static const bool debug_output = false;
  static const size_t max_motif_len = 31;
  static const size_t n_generalizations = 3 * max_motif_len + 1;

  matrix_t z(max_motif_len + 1, n_generalizations);
  for (size_t n = 0; n <= max_motif_len; n++) {
    for (size_t k = 0; k < n_generalizations; k++) {
      if (k == 0)
        z(n, k) = (1l << (2 * n));
      else if (3 * n < k)
        z(n, k) = 0;
      else {
        z(n, k) = 0;
        // to determine the number of IUPAC sequences of length n, we consider
        // the previously computed numbers of IUPAC sequences of length n-1 that
        //  A) have 3 degrees of degeneracy less than k
        //  B) have 2 degrees of degeneracy less than k
        //  C) have 1 degrees of degeneracy less than k
        //  D) have 0 degrees of degeneracy less than k
        // Each of these numbers of IUPAC sequences is multiplied with the
        // number of IUPAC symbols that, by being append to those sequences, add
        // the desired degree of degeneracy.
        // In the general case all four (A,B,C,D) are possible extensions.
        // However, when k is small then some of the added degrees of degeneracy
        // would exceed k.
        // Note that due to the absence of a break statement all following cases
        // are respected.
        switch (k) {
          default:                          // A
            z(n, k) += 1 * z(n - 1, k - 3); // add 3 degrees of degeneracy: one
                                            // IUPAC symbols allowing four
                                            // nucleotides: N
          case 2:                           // B
            z(n, k) += 4 * z(n - 1, k - 2); // add 2 degrees of degeneracy:
                                            // four IUPAC symbols allowing
                                            // three nucleotides: B, D, H, V
          case 1:                           // C
            z(n, k) += 6 * z(n - 1, k - 1); // add 1 degree of degeneracy: six
                                            // IUPAC symbols allowing two
                                            // nucleotides: W, S, M, K, R, Y
          case 0:                           // D
            z(n, k) += 4 * z(n - 1, k);     // add 0 degrees of degeneracy: four
                                            // IUPAC symbols allowing one
                                            // nucleotides: A, C, G, T
        }
      }
    }
  }

  if (debug_output) {
    cout << "Zeta table before accumulation" << endl;
    for (size_t n = 0; n <= max_motif_len; n++) {
      cout << n;
      for (size_t k = 0; k < n_generalizations; k++)
        cout << "\t" << z(n, k);
      cout << endl;
    }
  }

  // accumulate
  for (size_t n = 0; n <= max_motif_len; n++)
    for (size_t k = 1; k < n_generalizations; k++)
      z(n, k) += z(n, k - 1);

  if (debug_output) {
    cout << "Zeta table" << endl;
    for (size_t n = 0; n <= max_motif_len; n++) {
      cout << n;
      for (size_t k = 0; k < n_generalizations; k++)
        cout << "\t" << z(n, k);
      cout << endl;
    }
  }
  return z;
};

static const matrix_t zeta_table = compute_zeta_table();

double compute_correction(size_t length, size_t degeneracy) {
  double z = log(zeta_table(length, degeneracy));
  // cout << "log(zeta(" << length << ", " << degeneracy << ") = " << z << endl;
  return z;
}
