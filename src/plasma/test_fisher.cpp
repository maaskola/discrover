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

using namespace std;

int main(int argc, const char **argv)
{
  const size_t n = 50;

  matrix_t m(2,2);
  for(size_t i = 0; i < n; i++)
    for(size_t j = 0; j < n; j++)
      for(size_t k = 0; k < n; k++)
        for(size_t l = 0; l < n; l++) {
          m(0,0) = i;
          m(0,1) = j;
          m(1,0) = k;
          m(1,1) = l;
          cout << i << " " << j << " " << " " << k << " " << l << " " << fisher_exact_test(m).p_value << endl;
        }

  return EXIT_SUCCESS;
}

