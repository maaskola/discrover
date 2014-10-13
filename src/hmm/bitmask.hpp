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

#include <bitset>
#include <vector>

const static size_t max_motifs = 32;
typedef std::bitset<max_motifs> bitmask_t;

bitmask_t make_mask(const std::vector<size_t> &v);
std::vector<size_t> unpack_mask(const bitmask_t x);

#endif

