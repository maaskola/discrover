/*
 * =====================================================================================
 *
 *       Filename:  bitmask.hpp
 *
 *    Description:  Using bit sets to flag sets of motif indices
 *
 *        Created:  10/12/2014 08:34:55 AM
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef BITMASK_HPP
#define BITMASK_HPP

#include <stdexcept>
#include <bitset>
#include <vector>
#include <climits>

static const size_t max_motifs = sizeof(size_t) * CHAR_BIT;
using bitmask_t = std::bitset<max_motifs>;

bitmask_t make_mask(const std::vector<size_t> &v);
std::vector<size_t> unpack_mask(const bitmask_t x);

namespace Exception {
namespace BitMask {
struct TooManyMotifs : public std::runtime_error {
  TooManyMotifs(size_t num_motifs);
};
}
}

#endif
