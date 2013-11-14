/* =====================================================================================
 * Copyright (c) 2012, Jonas Maaskola
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
 *       Filename:  plasma_cli.hpp
 *
 *    Description:  Routine to construct the plasma CLI options
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef  PLASMA_CLI_HPP
#define  PLASMA_CLI_HPP

#include <boost/program_options.hpp>
#include "options.hpp"

boost::program_options::options_description gen_iupac_options_description(Seeding::options_t &options, const std::string &prefix="", const std::string &name="IUPAC regular expression finding options", size_t cols=80, bool include_all=true, bool allow_short=true);

#endif   /* ----- #ifndef PLASMA_CLI_HPP  ----- */

