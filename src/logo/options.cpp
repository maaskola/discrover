/*
 * =====================================================================================
 *
 *       Filename:  options.cpp
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
#include "options.hpp"

namespace Logo {
Options::Options()
    : pdf_logo(true),
      png_logo(true),
      axes(false),
      revcomp(false),
      type(Type::Sequence),
      alphabet(Alphabet::Undefined),
      order(Order::Frequency),
      palette(Palette::Default),
      absent(0.03),
      scale(1.0){}
}
