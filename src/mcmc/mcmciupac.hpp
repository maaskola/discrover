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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef MCMCHMM_HPP
#define MCMCHMM_HPP

#include "../plasma/score.hpp"
// #include "../plasma/motif.hpp"
#include "montecarlo.hpp"

namespace MCMC {
  typedef std::string Motif;
  template <>
    class Generator<Motif> {
      private:
        Seeding::Options options;
        size_t max_degeneracy;
//        size_t min_size, max_size;
        void replace_similar(char &c) const {
          if(options.verbosity >= Verbosity::debug)
            std::cout << "Replacing nucleotide " << static_cast<char>(c) << std::endl;
          switch(tolower(c)) {
            // individual nucleotides
            case 'a': c = "cgtmwr"[rand() % 6]; break;
            case 'c': c = "agtmsy"[rand() % 6]; break;
            case 'g': c = "actksr"[rand() % 6]; break;
            case 't':
            case 'u': c = "acgkwy"[rand() % 6]; break;
              // two nucleotide wildcards
            case 'm': c = "acrwsyvh"[rand() % 8]; break;
            case 'r': c = "agmwskvd"[rand() % 8]; break;
            case 'w': c = "atmrykhd"[rand() % 8]; break;
            case 's': c = "cgmyrkvb"[rand() % 8]; break;
            case 'y': c = "ctmswkhb"[rand() % 8]; break;
            case 'k': c = "gtmswydb"[rand() % 8]; break;
              // three nucleotide wildcards
            case 'b': c = "sykdhvn"[rand() % 7]; break;
            case 'd': c = "rwkbhvn"[rand() % 7]; break;
            case 'h': c = "mwykbvn"[rand() % 7]; break;
            case 'v': c = "mrsbdhn"[rand() % 7]; break;
              // four nucleotide wildcard
            case 'n': c = "bdhv"[rand() % 4]; break;
            default:
              throw("Error: character not recognized.");
          }
        }
        void replace_arbitrary(char &c) const {
          size_t r = rand() % 14;
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
        Generator(const Seeding::Options &opt, size_t w, size_t max_degen) : // , size_t min_size_=-1, size_t max_size_=-1) :
          options(opt), max_degeneracy(max_degen) // , min_size(min_size_), max_size(max_size_)
      {
        /*
        if(min_size == -1)
          min_size = w;
        if(max_size == -1)
          max_size = w;
  */
      };
        Motif generate(const Motif &motif_) const {
          const size_t R = 3;
          Motif motif(motif_);
          if(options.verbosity >= Verbosity::verbose)
            std::cout << "Generating new motif based off of " << motif << std::endl;
          do {
            size_t r = rand() % R;
            size_t p;
            switch(r) {
              case 0: // replace nucleotide by a similar one
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Replace nucleotide by a similar one" << std::endl;
                p = rand() % motif.size();
                replace_similar(motif[p]);
                break;
              case 1: // replace nucleotide by a random one
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Replace nucleotide by an arbitrary one" << std::endl;
                p = rand() % motif.size();
                replace_arbitrary(motif[p]);
                break;
              case 2: // roll one position
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "Roll one position" << std::endl;
                {
                  char nucl = "acgtmrwsykbdhvn"[rand() % 15];
                  std::string w = " ";
                  w[0] = nucl;
                  bool r = rand() % 2;
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
    };
  template <>
    class Evaluator<Motif> {
      private:
        Seeding::DataCollection data;
        Seeding::Options options;
        Seeding::Objective objective;
      public:
        Evaluator(const Seeding::DataCollection &data_, const Seeding::Options &opt, const Seeding::Objective &obj) :
          data(data_), options(opt), objective(obj) { };

        double evaluate(const Motif &motif) const {
          count_vector_t counts = count_motif(data, motif, options);
          double score = compute_score(data, counts, options, objective, motif.size(), Seeding::motif_degeneracy(motif));
          return(score);
        };
    };
}

 
#endif
