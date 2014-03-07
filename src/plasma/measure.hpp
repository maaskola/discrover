/* =====================================================================================
 * Copyright (c) 2012, Jonas Maaskola
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
 *       Filename:  measure.hpp
 *
 *    Description:  Discriminative objective functions
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <string>
#include <iostream>

namespace Measures {
  namespace Discrete {
    enum class Measure {
      Undefined,
      SignalFrequency,
      ControlFrequency,
      MutualInformation,
      VarianceMutualInformation,
      Gtest,
      LogpGtest,
      CorrectedLogpGtest,
      MatthewsCorrelationCoefficient,
      DeltaFrequency,
    };

    void parse_measure(const std::string &token, Measure &measure);
    std::string measure2string(Measure measure);
    std::string measure2acronym(Measure measure);

    std::istream &operator>>(std::istream &in, Measure &measure);
    std::ostream &operator<<(std::ostream &out, const Measure &measure);
  }

  namespace Continuous {
    enum class Measure {
      Likelihood,
      Viterbi,
      MutualInformation,
      ResidualMutualInformation,
      RankInformation,
      MatthewsCorrelationCoefficient,
      ChiSquare,
      LogLikelihoodDifference,
      DeltaFrequency,
      ClassificationPosterior,
      ClassificationLikelihood,
      Undefined
    };

    void parse_measure(const std::string &token, Measure &measure);
    std::string measure2string(Measure measure);

    std::istream& operator>>(std::istream& in, Measure &measure);
    std::ostream &operator<<(std::ostream &out, const Measure &measure);
  }

  template <typename X> bool is_discriminative(X measure);
  template <> bool is_discriminative<Discrete::Measure>(Discrete::Measure measure);
  template <> bool is_discriminative<Continuous::Measure>(Continuous::Measure measure);
  template <typename X> bool is_generative(X measure) { return(false); };
  template <> bool is_generative<Discrete::Measure>(Discrete::Measure measure);
  template <> bool is_generative<Continuous::Measure>(Continuous::Measure measure);
  template <typename X> bool is_two_by_two(X measure);
  template <> bool is_two_by_two<Discrete::Measure>(Discrete::Measure measure);
  template <> bool is_two_by_two<Continuous::Measure>(Continuous::Measure measure);
//  template <typename X> bool is_generative(X measure) { return(not is_discriminative<X>(measure)); }

  Discrete::Measure corresponding_measure(Continuous::Measure);
}

#endif   /* ----- #ifndef MEASURE_HPP  ----- */

