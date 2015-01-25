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
 *       Filename:  specification.hpp
 *
 *    Description:  Data structure to specify motifs that are to be searched
 *
 *        Created:  Fri Jun 22 02:08:38 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef SPECIFICATION_HPP
#define SPECIFICATION_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <set>

namespace Measures {
template <typename X>
bool is_discriminative(X measure);
template <typename X>
bool is_two_by_two(X measure);
};

namespace Specification {
namespace Exception {
struct MotifNamedControl : public std::runtime_error {
MotifNamedControl(const std::string &token);
};
struct ControlWithMotif : public std::runtime_error {
ControlWithMotif(const std::string &token);
};
}
/** the format is [MIDs[:CONTRAST:]]PATH
 * where
 * MIDs a set of motif ID (optional)
 * CONTRAST is the name of a contrast this data set belongs to; if unspecified,
 * the unnamed general contrast is used
 * PATH a path to a FASTA file
 */
struct Set {
  std::string contrast;
  std::string path;
  bool is_shuffle;
  bool is_control;
  std::set<std::string> motifs;
  Set(const std::string &s, bool shuffled = false);
  Set(const Set &spec);
  Set();
  std::string name() const;
};

using Sets = std::vector<Set>;

/** the format is [MID[:INSERT:]]MSPEC where
 * MID is a motif ID (optional)
 * INSERT is list of comma separated numbers that give positions after which insert states are allowed
 * MSPEC is one of three possibilities (tested in turn)
 *  1. a path to a file. If such a file exists, it will be opened and a PWM will be read from it
 *  2. a IUPAC regular expression
 *  3. a length specification for lib Plasma, optionally specifying multiplicity
 * examples:
 *   signal:2,3:5-9x3
 *     - the motif ID is "signal"
 *     - insertions after postions 2 and 3
 *     - lengths of between 5 and 9
 *     - three variants to consider
 *   previous:/home/user/motif.hmm
 *     - consider the file given above and name the motif "previous"
 *   motif:4:tgtanata
 *     - motif ID is "motif"
 *     - it is the IUPAC string tgtanata
 *     - it allows one insertion after the first "a"
 */
struct Motif {
  enum class Kind { File, Seed, Plasma };
  Kind kind;
  std::string specification;
  std::string name;
  std::vector<size_t> insertions;
  std::vector<size_t> lengths;
  size_t multiplicity;
  Motif(const std::string &s = "");
  Motif(const Motif &spec);
};

using Motifs = std::vector<Motif>;

namespace Contrast {
/**
 * CONTRASTEXP is a contrast expression, as below
 * CONTRASTEXP: ATOM | CONTRASTEXP-ATOM | CONTRASTEXP+ATOM
 */
struct Atom {
  int sign;
  std::string contrast;
};

using Expression = std::vector<Atom>;
Expression make_contrast(const std::string &path);

std::string to_string(const Atom &atom);
std::string to_string(const Expression &expr);
std::ostream &operator<<(std::ostream &out, const Atom &atom);
std::ostream &operator<<(std::ostream &out, const Expression &expr);
}

/** the format is [MID[:CONTRASTEXP:]]OBJ where
 * MID is a motif ID (optional)
 *  CONTRASTEXP is a contrast expression, as below
 *  OBJ is a string parsable as the desired objection function type (of template type X),
 *  i.e. it requires a definition of void parse_measure(const std::string&, &X)
 *  TODO: document that contrast names may not contain - or +
 */
template <typename X>
struct Objective {
  using measure_t = X;
  using contrast_atom_t = Contrast::Atom;
  using contrast_expression_t = Contrast::Expression;
  using iterator = typename contrast_expression_t::iterator;
  using const_iterator = typename contrast_expression_t::const_iterator;

  std::string motif_name;
  contrast_expression_t contrast_expression;
  measure_t measure;

  Objective()
      : motif_name(""), contrast_expression(), measure(measure_t::Undefined){};

  Objective(const std::string &s)
      : motif_name(""), contrast_expression(), measure() {
    std::string specification = s;
    size_t pos;
    if ((pos = specification.find(":")) != std::string::npos) {
      motif_name = specification.substr(0, pos);
      specification = specification.substr(pos + 1);
    }
    if ((pos = specification.find(":")) != std::string::npos) {
      std::string expr_s = specification.substr(0, pos);
      contrast_expression = Contrast::make_contrast(expr_s);
      specification = specification.substr(pos + 1);
    }
    parse_measure(specification, measure);
  };

  Objective(const Objective &obj)
      : motif_name(obj.motif_name),
        contrast_expression(obj.contrast_expression),
        measure(obj.measure){};
};

template <typename X>
typename Objective<X>::iterator begin(Objective<X> &objective) {
  return begin(objective.contrast_expression);
}
template <typename X>
typename Objective<X>::iterator end(Objective<X> &objective) {
  return end(objective.contrast_expression);
}
template <typename X>
typename Objective<X>::const_iterator begin(const Objective<X> &objective) {
  return begin(objective.contrast_expression);
}
template <typename X>
typename Objective<X>::const_iterator end(const Objective<X> &objective) {
  return end(objective.contrast_expression);
}

std::istream &operator>>(std::istream &in, Motif &m);
std::istream &operator>>(std::istream &in, Set &spec);
template <typename X>
std::istream &operator>>(std::istream &in, Objective<X> &spec) {
  std::string token;
  in >> token;
  spec = Objective<X>(token);
  return in;
}

std::string to_string(const Motif &spec);
std::string to_string(const Set &spec);
template <typename X>
std::string to_string(const Objective<X> &spec) {
  std::string s = std::string(spec.motif_name) + ":"
                  + to_string(spec.contrast_expression) + ":"
                  + measure2string(spec.measure);
  return s;
}

std::ostream &operator<<(std::ostream &out, const Motif &m);
std::ostream &operator<<(std::ostream &out, const Set &spec);
template <typename X>
std::ostream &operator<<(std::ostream &out, Objective<X> &spec) {
  out << to_string(spec);
  return out;
}
}
#endif /* ----- #ifndef SPECIFICATION_HPP ----- */
