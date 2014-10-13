/*
 * =====================================================================================
 *
 *       Filename:  results.cpp
 *
 *    Description:  
 *
 *        Created:  Thu Mar 07 20:30:35 2014 +0200
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "results.hpp"

namespace Training {
  State::State(size_t n) : center(-9), scores(n) { };

  Result::Result() : state(), delta(0), parameter_file("") { };
}

