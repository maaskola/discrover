/*
 * =====================================================================================
 *
 *       Filename:  measure.hpp
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

/*namespace Plasma {
  // typedef Measures::Discrete::Measure Measure;
  namespace HMM {
    typedef Measures::Continuous::Measure Measure;
  }
}*/

#endif   /* ----- #ifndef MEASURE_HPP  ----- */

