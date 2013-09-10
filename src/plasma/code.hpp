/*
 * =====================================================================================
 *
 *       Filename:  code.hpp
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


#ifndef  CODE_HPP
#define  CODE_HPP

#include <string>

namespace Plasma {
  const char* const Symbol = "-acmgrsvtwyhkdbn";

  char* construct_code();

  extern char* Code;

  std::string iupac2regex(const std::string &s);
};

#endif   /* ----- #ifndef CODE_HPP ----- */

