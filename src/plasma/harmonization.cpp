/*
 * =====================================================================================
 *
 *       Filename:  harmonization.cpp
 *
 *    Description:  Ensure that motif, data, and objective specifications are consistent
 *
 *        Created:  Sun Jan 25 2015 16:05:34 +0200
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "harmonization.hpp"

using namespace std;

namespace Specification {
namespace Exception {
namespace Motif {
NameNotUnique::NameNotUnique(const std::string &name)
    : runtime_error("Error: motif name '" + name + "' not unique.") {}
MultiplicityZero::MultiplicityZero(const std::string &name)
    : runtime_error("Error: multiplicity of motif '" + name + "' is zero!") {}
NameNotUniqueInObjective::NameNotUniqueInObjective(const std::string &name)
    : runtime_error("Error: motif name '" + name
                    + "'in objective not unique.") {}
NoSpecfication::NoSpecfication(const std::string &name)
    : runtime_error(
          "Error: found a motif name in the objective for which "
          "no motif specification is found: '" + name + "'.") {}
WhenOneThenAll::WhenOneThenAll()
    : runtime_error(
          "Error: when any objectives name motifs then all "
          "objectives must define motifs.") {}
}
namespace Objective {
MultpleObjectivesWithoutNamedMotifs::MultpleObjectivesWithoutNamedMotifs()
    : runtime_error(
          "Error: multiple score specifications given that do not name "
          "motifs.") {}
}
}
}
