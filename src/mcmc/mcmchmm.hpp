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
 *       Filename:  mcmchmm.hpp
 *
 *    Description:  Code for MCMC sampling of HMM parameters
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef MCMCHMM_HPP
#define MCMCHMM_HPP

#include "../hmm/hmm.hpp"
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
          return(hmm.random_variant(options, EntropySource::rng));
        };
    };
  template <>
    class Evaluator<HMM> {
      private:
        Data::Collection data;
        Training::Task task;
        Training::State ts;
        hmm_options options;
      public:
        Evaluator(const Data::Collection &data_,
            const Training::Task &obj,
            const hmm_options &opt) :
          data(data_),
          task(obj),
          ts(1), // TODO make this work with multiple tasks
          options(opt) {
          };

        double evaluate(const HMM &hmm) const {
          return(hmm.compute_score(data, task, options.weighting));
        };
    };
}

 
#endif
