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
 *       Filename:  mcmciupac.hpp
 *
 *    Description:  Code for MCMC sampling of IUPAC regular expression motifs
 *
 *        Created:  Wed Nov 13 21:28:55 2013 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef MCMCHMM_HPP
#define MCMCHMM_HPP

#include <random>
#include <functional>
#include "../plasma/score.hpp"
#include "../plasma/motif.hpp"
#include "montecarlo.hpp"

namespace MCMC {
  using Motif = std::string;
  template <>
    class Generator<Motif> {
      private:
        Seeding::Options options;
        size_t motif_length;
        size_t max_degeneracy;
        std::uniform_int_distribution<size_t> dist2, dist3, dist4, dist6, dist7, dist8, dist14, dist15, distPos;
        std::function<size_t()>               rand2, rand3, rand4, rand6, rand7, rand8, rand14, rand15, randPos;
        void replace_similar(char &c) const {
          if(options.verbosity >= Verbosity::debug)
            std::cout << "Replacing nucleotide " << static_cast<char>(c) << std::endl;
          switch(tolower(c)) {
            // individual nucleotides
            case 'a': c = "cgtmwr"[rand6()]; break;
            case 'c': c = "agtmsy"[rand6()]; break;
            case 'g': c = "actksr"[rand6()]; break;
            case 't':
            case 'u': c = "acgkwy"[rand6()]; break;
              // two nucleotide wildcards
            case 'm': c = "acrwsyvh"[rand8()]; break;
            case 'r': c = "agmwskvd"[rand8()]; break;
            case 'w': c = "atmrykhd"[rand8()]; break;
            case 's': c = "cgmyrkvb"[rand8()]; break;
            case 'y': c = "ctmswkhb"[rand8()]; break;
            case 'k': c = "gtmswydb"[rand8()]; break;
              // three nucleotide wildcards
            case 'b': c = "sykdhvn"[rand7()]; break;
            case 'd': c = "rwkbhvn"[rand7()]; break;
            case 'h': c = "mwykbvn"[rand7()]; break;
            case 'v': c = "mrsbdhn"[rand7()]; break;
              // four nucleotide wildcard
            case 'n': c = "bdhv"[rand4()]; break;
            default:
              throw("Error: character not recognized.");
          }
        }
        void replace_arbitrary(char &c) const {
          const size_t r = rand14();
          switch(tolower(c)) {
            // individual nucleotides
            case 'a': c = "cgtmrwsykbdhvn"[r]; break;
            case 'c': c = "agtmrwsykbdhvn"[r]; break;
            case 'g': c = "actmrwsykbdhvn"[r]; break;
            case 't':
            case 'u': c = "acgmrwsykbdhvn"[r]; break;
            case 'm': c = "acgtrwsykbdhvn"[r]; break;
            case 'r': c = "acgtmwsykbdhvn"[r]; break;
            case 'w': c = "acgtmrsykbdhvn"[r]; break;
            case 's': c = "acgtmrwykbdhvn"[r]; break;
            case 'y': c = "acgtmrwskbdhvn"[r]; break;
            case 'k': c = "acgtmrwsybdhvn"[r]; break;
            case 'b': c = "acgtmrwsykdhvn"[r]; break;
            case 'd': c = "acgtmrwsykbhvn"[r]; break;
            case 'h': c = "acgtmrwsykbdvn"[r]; break;
            case 'v': c = "acgtmrwsykbdhn"[r]; break;
            case 'n': c = "acgtmrwsykbdhv"[r]; break;
            default: throw("Error: character not recognized.");
          }
        }
      public:
        Generator(const Seeding::Options &opt, size_t len, size_t max_degen) :
          options(opt),
          motif_length(len),
          max_degeneracy(max_degen),
          dist2(0, 1), dist3(0, 2), dist4(0, 3), dist6(0, 5),
          dist7(0, 6), dist8(0, 7), dist14(0, 13), dist15(0, 14),
          distPos(0, len-1) {
            rand2 = std::bind(dist2, MCMC::EntropySource::rng);
            rand3 = std::bind(dist3, MCMC::EntropySource::rng);
            rand4 = std::bind(dist4, MCMC::EntropySource::rng);
            rand6 = std::bind(dist6, MCMC::EntropySource::rng);
            rand7 = std::bind(dist7, MCMC::EntropySource::rng);
            rand8 = std::bind(dist8, MCMC::EntropySource::rng);
            rand14 = std::bind(dist14, MCMC::EntropySource::rng);
            rand15 = std::bind(dist15, MCMC::EntropySource::rng);
            randPos = std::bind(distPos, MCMC::EntropySource::rng);
          };
        Motif generate(const Motif &motif_) const {
          Motif motif(motif_);
          if(options.verbosity >= Verbosity::verbose)
            std::cout << "Generating new motif based off of " << motif << std::endl;
          do {
            size_t r = rand3();
            size_t p;
            switch(r) {
              case 0: // replace nucleotide by a similar one
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Replace nucleotide by a similar one" << std::endl;
                p = randPos();
                replace_similar(motif[p]);
                break;
              case 1: // replace nucleotide by a random one
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Replace nucleotide by an arbitrary one" << std::endl;
                p = randPos();
                replace_arbitrary(motif[p]);
                break;
              case 2: // roll one position
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Roll one position" << std::endl;
                {
                  char nucl = "acgtmrwsykbdhvn"[rand15()];
                  std::string w = " ";
                  w[0] = nucl;
                  bool r = rand2();
                  size_t n = motif.size() - 1;
                  if(r == 0)
                    motif = w + motif.substr(0, n);
                  else
                    motif = motif.substr(1, n) + w;
                }
                break;
            }
          } while(Seeding::motif_degeneracy(motif) > max_degeneracy);
          if(options.verbosity >= Verbosity::verbose)
            std::cout << motif_ << " -> " << motif << std::endl;
          return(motif);
        };
        Motif generate() {
          std::string word;
          for(size_t j = 0; j < motif_length; j++)
            word += "acgt"[rand4()];
          return(word);
        };
    };
  template <>
    class Evaluator<Motif> {
      private:
        Seeding::Collection collection;
        Seeding::Options options;
        Seeding::Objective objective;
      public:
        Evaluator(const Seeding::Collection &col, const Seeding::Options &opt, const Seeding::Objective &obj) :
          collection(col), options(opt), objective(obj) { };

        double evaluate(const Motif &motif) const {
          count_vector_t counts = count_motif(collection, motif, options);
          double score = compute_score(collection, counts, options, objective, motif.size(), Seeding::motif_degeneracy(motif));
          return(score);
        };
    };
}

#endif
