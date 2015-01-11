/*
 * =====================================================================================
 *
 *       Filename:  bitmask.cpp
 *
 *    Description:
 *
 *        Created:  10/12/2014 01:57:44 AM
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include <sstream>
#include "bitmask.hpp"

using namespace std;

bitmask_t make_mask(const vector<size_t> &v) {
  bitmask_t x = 0;
  for (auto y : v) {
    if (y >= max_motifs) {
      throw Exception::BitMask::TooManyMotifs(y);
    }
    // TODO implement in terms of set() or operator[] methods
    x |= 1 << y;
  }
  return x;
}

vector<size_t> unpack_mask(const bitmask_t x) {
  vector<size_t> v;
  size_t z = 0;
  bitmask_t y = 1;
  while (x.to_ullong() >= y.to_ullong()) {
    if ((x & y) != 0)
      v.push_back(z);
    // TODO implement in terms of test() or operator[] methods
    y = y << 1;
    z++;
  }
  return v;
}

namespace Exception {
namespace BitMask {
TooManyMotifs::TooManyMotifs(size_t num_motifs_)
    : exception(), num_motifs(num_motifs_) {};
const char *TooManyMotifs::what() const noexcept {
  stringstream ss;
  ss << "Error: trying to construct mask for too many motifs! "
     << "The offending index is: " << num_motifs
     << ", and this version only supports " << max_motifs << " motifs." << endl;
  return ss.str().c_str();
}
}
}
