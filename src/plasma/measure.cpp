/*
 * =====================================================================================
 *
 *       Filename:  measure.cpp
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

#include <iostream>
#include <cassert>
#include "measure.hpp"

using namespace std;

namespace Measures {
  namespace Discrete {
    void parse_measure(const string &token, Measure &measure) {
      if(token == "mi" or token == "mico")
        measure = Measure::MutualInformation;
      else if(token == "signal_freq")
        measure = Measure::SignalFrequency;
      else if(token == "control_freq")
        measure = Measure::ControlFrequency;
      else if(token == "mcc")
        measure = Measure::MatthewsCorrelationCoefficient;
      else if(token == "delta_freq" or token == "dfreq")
        measure = Measure::DeltaFrequency;
      else if(token == "gtest")
        measure = Measure::Gtest;
      else if(token == "gtest_logp_raw")
        measure = Measure::LogpGtest;
      else if(token == "gtest_logp")
        measure = Measure::CorrectedLogpGtest;
      else if(token == "fisher")
        measure = Measure::FisherExactTest;
      else
        measure = Measure::Undefined;
    }

    string measure2string(Measure measure) {
      switch(measure) {
        case Measure::SignalFrequency:
          return("signal relative frequency");
        case Measure::ControlFrequency:
          return("control relative frequency");
        case Measure::MutualInformation:
          return("mutual information");
        case Measure::MatthewsCorrelationCoefficient:
          return("Matthew's correlation coefficient");
        case Measure::DeltaFrequency:
          return("delta frequency");
        case Measure::Gtest:
          return("G-test");
        case Measure::LogpGtest:
          return("Uncorrected log P(G-test)");
        case Measure::CorrectedLogpGtest:
          return("Bonferroni-corrected log P(G-test)");
        case Measure::VarianceMutualInformation:
          return("variance of mutual information");
        case Measure::FisherExactTest:
          return("Fisher's exact test");
        case Measure::Undefined:
          return("undefined");
      }
      return("");
    }


    string measure2acronym(Measure measure) {
      switch(measure) {
        case Measure::SignalFrequency:
          return("signal_freq");
        case Measure::ControlFrequency:
          return("control_freq");
        case Measure::MutualInformation:
          return("mi");
        case Measure::MatthewsCorrelationCoefficient:
          return("mcc");
        case Measure::DeltaFrequency:
          return("delta_freq");
        case Measure::Gtest:
          return("gtest");
        case Measure::LogpGtest:
          return("gtest_logp_raw");
        case Measure::CorrectedLogpGtest:
          return("gtest_logp");
        case Measure::VarianceMutualInformation:
          return("var_mi");
        case Measure::FisherExactTest:
          return("fisher");
        case Measure::Undefined:
          return("undefined");
      }
      return("");
    };

    istream &operator>>(istream &in, Measure &measure) {
      string token;
      in >> token;
      parse_measure(token, measure);
      return(in);
    }

    ostream &operator<<(ostream &out, const Measure &measure) {
      out << measure2string(measure);
      return(out);
    }
  }

  template <> bool is_discriminative<Continuous::Measure>(Continuous::Measure measure)
  {
    if(measure == Continuous::Measure::Likelihood or measure == Continuous::Measure::Viterbi)
      return false;
    else
      return true;
  }

  template <> bool is_discriminative<Discrete::Measure>(Discrete::Measure measure)
  {
    if(measure == Discrete::Measure::Undefined
        or measure == Discrete::Measure::SignalFrequency
        or measure == Discrete::Measure::ControlFrequency)
      return false;
    else
      return true;
  }

  template <> bool is_generative<Continuous::Measure>(Continuous::Measure measure)
  {
    return(not is_discriminative(measure));
  };

  template <> bool is_two_by_two<Continuous::Measure>(Continuous::Measure measure)
  {
    if(measure == Continuous::Measure::DeltaFrequency
        or measure == Continuous::Measure::LogLikelihoodDifference
        or measure == Continuous::Measure::MatthewsCorrelationCoefficient)
      return true;
    else
      return false;
  }

  template <> bool is_two_by_two<Discrete::Measure>(Discrete::Measure measure)
  {
    if(measure == Discrete::Measure::MatthewsCorrelationCoefficient
        or measure == Discrete::Measure::DeltaFrequency)
      return true;
    else
      return false;
  }


  namespace Continuous {
    void parse_measure(const string &token, Measure &measure) {
      if(token == "likelihood")
        measure = Measure::Likelihood;
      else if(token == "viterbi")
        measure = Measure::Viterbi;
      else if(token == "bw")
        measure = Measure::Likelihood;
      else if(token == "mi" or token == "mico")
        measure = Measure::MutualInformation;
      else if(token == "ri")
        measure = Measure::RankInformation;
      else if(token == "mcc")
        measure = Measure::MatthewsCorrelationCoefficient;
      else if(token == "chisq")
        measure = Measure::ChiSquare;
      else if(token == "dlogl")
        measure = Measure::LogLikelihoodDifference;
      else if(token == "dme") {
        cout << "Deprecation warning: 'dme' has been renamed to 'dlogl'. Consider using that instead, as support for 'dme' will be removed in the future." << endl;
        measure = Measure::LogLikelihoodDifference;
      } else if(token == "class" or token == "mmie")
        measure = Measure::ClassificationPosterior;
      else if(token == "class-like" or token == "mmie-like")
        measure = Measure::ClassificationLikelihood;
      else if(token == "dfreq")
        measure = Measure::DeltaFrequency;
      else if(token == "dips") {
        cout << "Deprecation warning: 'dips' has been renamed to 'dfreq'. Consider using that instead, as support for 'dips' will be removed in the future." << endl;
        measure = Measure::DeltaFrequency;
      } else if(token == "none" or token == "undefined") {
        measure = Measure::Undefined;
      } else {
        cout << "Measure '" << token << "' not implemented. See -h or --help for help." << endl;
        assert(false);
        exit(-1);
      }
    }

    istream& operator>>(istream& in, Measure& measure)
    {
      string token;
      in >> token;
      parse_measure(token, measure);
      return(in);
    }

    ostream &operator<<(ostream &out, const Measure &measure) {
      out << measure2string(measure);
      return(out);
    }

    string measure2string(Measure measure)
    {
      string s;
      switch(measure) {
        case Measure::Likelihood:
          s = "likelihood";
          break;
        case Measure::Viterbi:
          s = "viterbi";
          break;
        case Measure::MutualInformation:
          s = "mi";
          break;
        case Measure::RankInformation:
          s = "ri";
          break;
        case Measure::MatthewsCorrelationCoefficient:
          s = "mcc";
          break;
        case Measure::ChiSquare:
          s = "chisq";
          break;
        case Measure::LogLikelihoodDifference:
          s = "dlogl";
          break;
        case Measure::ClassificationPosterior:
          // s = "class";
          s = "mmie";
          break;
        case Measure::ClassificationLikelihood:
          // s = "class-like";
          s = "mmie-like";
          break;
        case Measure::DeltaFrequency:
          s = "dfreq";
          break;
        case Measure::Undefined:
          s = "undefined";
          break;
      }
      return(s);
    }

  }

  Discrete::Measure corresponding_measure(Continuous::Measure hmm_measure) {
    switch(hmm_measure) {
      case Continuous::Measure::Likelihood:
      case Continuous::Measure::Viterbi:
        return(Discrete::Measure::SignalFrequency);
      case Continuous::Measure::MutualInformation:
      case Continuous::Measure::RankInformation:
      case Continuous::Measure::LogLikelihoodDifference:
      case Continuous::Measure::ClassificationPosterior:
      case Continuous::Measure::ClassificationLikelihood:
        return(Discrete::Measure::MutualInformation);
      case Continuous::Measure::MatthewsCorrelationCoefficient:
        return(Discrete::Measure::MatthewsCorrelationCoefficient);
      case Continuous::Measure::DeltaFrequency:
        return(Discrete::Measure::DeltaFrequency);
      case Continuous::Measure::ChiSquare:
      case Continuous::Measure::Undefined:
        return(Discrete::Measure::Undefined);
      default:
        cout << "Error: can't determine a corresponding measure for " << hmm_measure << endl;
        exit(-1);
    }
  }
}

