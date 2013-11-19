/*
 * =====================================================================================
 *
 *       Filename:  test_fisher.cpp
 *
 *    Description:  Tests the implementation of Fisher's exact test
 *
 *        Version:  1.0
 *        Created:  18.11.2013 19:55:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include "plasma_stats.hpp"
#include "dreme_hyper.hpp"

using namespace std;

int main(int argc, const char **argv)
{
  const size_t n = 50;

  const double log_pthresh = 0;
  matrix_t m(2,2);
  for(size_t i = 10; i < 10 + n; i++)
    for(size_t j = 1000; j < 1000 + n; j++)
      for(size_t k = 0; k < n; k++)
        for(size_t l = 10000; l < 10000 + n; l++) {
          m(0,0) = i;
          m(0,1) = j;
          m(1,0) = k;
          m(1,1) = l;
          auto x = fisher_exact_test(m);
          auto y = fisher_exact_test(m, 1, Alternative::Greater);
          double dreme_logp = getLogFETPvalue(i, i+j, k, k+l, log_pthresh);
          cout << i
            << "\t" << j
            << "\t" << k
            << "\t" << l
            << "\t" << x.p_value
            << "\t" << x.log_p_value
            << "\t" << y.p_value
            << "\t" << y.log_p_value
            << "\t" << dreme_logp
            << endl;
        }

  return EXIT_SUCCESS;
}

