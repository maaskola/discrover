/*
 * =====================================================================================
 *
 *       Filename:  results.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */

#include "results.hpp"

namespace Seeding {
  Result::Result(const Objective &objective) : Objective(objective), motif(), score(), log_p(), counts() { };
}

