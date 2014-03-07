/*
 * =====================================================================================
 *
 *       Filename:  results.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Thu Mar 07 20:30:35 2014 +0200
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */

#include "results.hpp"

namespace Training {
  State::State(size_t n) : center(-9), scores(n) { };

  Result::Result() : state(), delta(0), parameter_file("") { };
}

