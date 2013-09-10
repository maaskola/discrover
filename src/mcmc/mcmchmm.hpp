
#ifndef MCMCHMM_HPP
#define MCMCHMM_HPP

#include "../hmm.hpp"
#include "montecarlo.hpp"

namespace MCMC {
  template <>
    class Generator<HMM> {
      private:
        hmm_options options;
      public:
        Generator(const hmm_options &opt, size_t w) :
          options(opt)
        {
          if(options.sampling.min_size == -1)
            options.sampling.min_size = w;
          if(options.sampling.max_size == -1)
            options.sampling.max_size = w;
        };
        HMM generate(const HMM &hmm) const {
          return(hmm.random_variant(options));
        };
    };
  template <>
    class Evaluator<HMM> {
      private:
        Data::Collection data;
        Training::Task task;
        Training::State ts;
      public:
        Evaluator(const Data::Collection &data_,
            const Training::Task &obj) :
          data(data_),
          task(obj),
          ts(1) { // TODO make this work with multiple tasks
          };

        double evaluate(const HMM &hmm) const {
          return(hmm.compute_score(data, task));
        };
    };
}

 
#endif
