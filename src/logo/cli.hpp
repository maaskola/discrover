/* =====================================================================================
 * Copyright (c) 2014, Jonas Maaskola
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
 *       Filename:  cli.hpp
 *
 *    Description:  Routine to construct the CLI options for logo creation
 *
 *        Created:  Fri Jan 02 20:32:15 2014 +0100
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef LOGO_CLI_HPP
#define LOGO_CLI_HPP

#include <boost/program_options.hpp>
#include "options.hpp"

boost::program_options::options_description gen_logo_options_description(
    Logo::Options &options, bool iupac_mode = true, size_t cols = 80,
    const std::string &name = "Sequence logo creation options");

#endif /* ----- #ifndef LOGO_CLI_HPP  ----- */
