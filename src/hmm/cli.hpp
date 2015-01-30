/* =====================================================================================
 * Copyright (c) 2015, Jonas Maaskola
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
 *    Description:  Command line options generation for Discrover
 *
 *        Created:  Fri Jan 30 20:40:22 2015 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */
#ifndef DISCROVER_CLI_HPP
#define DISCROVER_CLI_HPP
#include <string>
#include <boost/program_options.hpp>
#include "hmm_options.hpp"

void gen_discrover_cli(
    size_t cols, std::string &config_path, Options::HMM &options,
    boost::program_options::options_description &common_options,
    boost::program_options::options_description &visible_options,
    boost::program_options::options_description &hidden_options,
    boost::program_options::options_description &cmdline_options,
    boost::program_options::options_description &config_file_options);
#endif
