/*
 * =====================================================================================
 *
 *       Filename:  bitmask.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/12/2014 01:57:44 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola (jonas@maaskola.de)
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>
#include "bitmask.hpp"

using namespace std;

bitmask_t make_mask(const vector<size_t> &v) {
  bitmask_t x = 0;
  for(auto y: v) {
    if(y > max_motifs) {
      cout << "Error: trying to construct mask for too many motifs! The offending index was: " << y << ", and there may only be " << max_motifs << " motifs in this version." << endl;
      exit(-1);
    }
    // TODO implement in terms of set() or operator[] methods
    x |= 1 << y;
  }
  return(x);
}

vector<size_t> unpack_mask(const bitmask_t x) {
  vector<size_t> v;
  size_t z = 0;
  bitmask_t y = 1;
  while(x.to_ullong() >= y.to_ullong()) {
    if((x & y) != 0)
      v.push_back(z);
    // TODO implement in terms of test() or operator[] methods
    y = y << 1;
    z++;
  }
  return(v);
}

