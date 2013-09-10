/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
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

#ifndef  PLASMA_CLI_HPP
#define  PLASMA_CLI_HPP

#include <boost/program_options.hpp>
#include "options.hpp"

boost::program_options::options_description gen_iupac_options_description(Plasma::options_t &options, const std::string &prefix="", const std::string &name="IUPAC regular expression finding options", size_t cols=80, bool include_all=true, bool allow_short=true);

#endif   /* ----- #ifndef PLASMA_CLI_HPP  ----- */

