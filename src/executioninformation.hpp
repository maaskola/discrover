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
 *       Filename:  executioninformation.hpp
 *
 *    Description:  Structure storing information about the current program execution
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef EXECUTIONINFORMATION_HPP
#define EXECUTIONINFORMATION_HPP

#include <string>

struct ExecutionInformation {
  std::string program_name;
  std::string hmm_version;
  std::string cmdline;
  std::string datetime;
  std::string directory;
};

std::string cmdline(int argc, const char** argv);
ExecutionInformation generate_exec_info(const std::string &name, const std::string &hmm_version, const std::string &cmdline);

#endif
