/*
 * =====================================================================================
 *
 *       Filename:  options.hpp
 *
 *    Description:  Option for sequence logo plotting
 *
 *        Version:  1.0
 *        Created:  01/02/2015 08:36:47 PM
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef LOGO_OPTIONS_HPP
#define LOGO_OPTIONS_HPP

namespace Logo {
enum class Type { Sequence, Frequency };

enum class Alphabet { RNA, DNA, Undefined };

enum class Order { Alphabetic, Frequency };

enum class Palette { Default, Tetrad, Solarized };

struct Options {
  Options();
  bool pdf_logo;
  bool png_logo;
  bool axes;
  bool revcomp;
  Type type;
  Alphabet alphabet;
  Order order;
  Palette palette;
  double absent;
  double scale;
};
}

#endif
