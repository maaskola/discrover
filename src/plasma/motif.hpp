/*
 * =====================================================================================
 *
 *       Filename:  motif.hpp
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


#ifndef  MOTIF_HPP
#define  MOTIF_HPP

#include <cstdint>
#include <list>
#include <string>
#include "data.hpp"
#include "stats.hpp"
#include "../verbosity.hpp"

namespace Plasma {
  /** A class for IUPAC regular expression type motifs of up to length 16.
   * Uses a binary representation of size 64 bit with 4 bits per position. */
/*  struct Motif {
    // TODO write long motif class
    // TODO use SSE2
    typedef std::string storage_t;
//    static const storage_t block = 0;
    storage_t data;
    size_t length;
    Motif();
    Motif(storage_t d);
    Motif(const std::string &s);

    Motif reverse_complement() const;
    std::string to_string(storage_t x) const;
    std::string to_string() const;

    bool match(const Sequence &seq, size_t &pos, Verbosity verbosity) const;

    size_t degeneracy() const;
  };
  */
  double information_content(const std::string &motif);
  size_t motif_degeneracy(const std::string &motif);
  template <typename Func> void each_generalization(const std::string &motif, Func func) {
    // std::cout << "Generalizations of " << to_string() << std::endl;
    std::string generalization = motif;
    for(size_t i = 0; i < motif.size(); i++) {
      char current = generalization[i];
      switch(current) {
        case 'a':
          generalization[i] = 'm';
          func(generalization);
          generalization[i] = 'r';
          func(generalization);
          generalization[i] = 'w';
          func(generalization);
          break;
        case 'c':
          generalization[i] = 'm';
          func(generalization);
          generalization[i] = 's';
          func(generalization);
          generalization[i] = 'y';
          func(generalization);
          break;
        case 'g':
          generalization[i] = 'r';
          func(generalization);
          generalization[i] = 's';
          func(generalization);
          generalization[i] = 'k';
          func(generalization);
          break;
        case 't':
          generalization[i] = 'w';
          func(generalization);
          generalization[i] = 'y';
          func(generalization);
          generalization[i] = 'k';
          func(generalization);
          break;
        case 'w':
          generalization[i] = 'd';
          func(generalization);
          generalization[i] = 'h';
          func(generalization);
        case 's':
          generalization[i] = 'b';
          func(generalization);
          generalization[i] = 'v';
          func(generalization);
        case 'm':
          generalization[i] = 'h';
          func(generalization);
          generalization[i] = 'v';
          func(generalization);
        case 'k':
          generalization[i] = 'b';
          func(generalization);
          generalization[i] = 'd';
          func(generalization);
        case 'r':
          generalization[i] = 'd';
          func(generalization);
          generalization[i] = 'v';
          func(generalization);
        case 'y':
          generalization[i] = 'b';
          func(generalization);
          generalization[i] = 'h';
          func(generalization);
        case 'b':
        case 'd':
        case 'h':
        case 'v':
          generalization[i] = 'n';
          func(generalization);
        default:
          break;
      }
      generalization[i] = current;
    }
  };

  std::list<std::string> all_generalizations(const std::string &motif);
};

// std::ostream &operator<<(std::ostream &os, const Plasma::Motif &motif);

#endif   /* ----- #ifndef MOTIF_HPP ----- */

