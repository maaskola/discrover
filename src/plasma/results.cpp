/*
 * =====================================================================================
 *
 *       Filename:  results.cpp
 *
 *    Description:  
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "results.hpp"

namespace Seeding {
  Result::Result(const Objective &objective) : Objective(objective), motif(), score(), log_p(), counts() { };
}

