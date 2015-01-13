/* =====================================================================================
 * Copyright (c) 2011, Jonas Maaskola
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 *
 *       Filename:  trainingmode.hpp
 *
 *    Description:  Types to represent learning method and objective function
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef TRAININGMODE_HPP
#define TRAININGMODE_HPP

#include <iostream>
#include <vector>
#include "../plasma/specification.hpp"
#include "../plasma/measure.hpp"
#include "../plasma/options.hpp"

using Measure = Measures::Continuous::Measure;

namespace Training {
enum class Method { None, Reestimation, Gradient };
std::string method2string(Method method);
std::istream& operator>>(std::istream& in, Method& method);
std::ostream& operator<<(std::ostream& out, Method method);

Method measure2method(Measure measure);

/** A Range represents a set of state indices */
using Range = std::vector<size_t>;

/** Targets represent a set of indices of states whose transition or emission
 * probabilities are to respected */
struct Targets {
  Range transition;
  Range emission;
};

using Objective = Specification::Objective<Measure>;
using Objectives = std::vector<Objective>;

Seeding::Objective corresponding_objective(const Objective& x,
                                           bool use_mi_to_seed);
Seeding::Objectives corresponding_objectives(const Objectives& x,
                                             bool use_mi_to_seed);

struct Task : public Objective {
  Targets targets;
};

using Tasks = std::vector<Task>;
};

namespace Exception {
namespace HMM {
struct InvalidTrainingMethod : public std::runtime_error {
  InvalidTrainingMethod(const std::string& token);
};
}
}

#endif
