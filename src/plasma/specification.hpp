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

#include <iostream> // TODO remove
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>

namespace Measures {
  template <typename X> bool is_discriminative(X measure);
  template <typename X> bool is_two_by_two(X measure);
};

namespace Specification {
  /** the format is [MIDs[:CONTRAST:]]PATH
   * where
   * MIDs a set of motif ID (optional)
   * CONTRAST is the name of a contrast this data set belongs to; if unspecified, the unnamed general contrast is used
   * PATH a path to a FASTA file
   */
  struct Set {
    std::string contrast;
    std::string path;
    bool is_shuffle;
    std::set<std::string> motifs;
    Set(const std::string &s, bool shuffled=false);
    Set(const Set &spec);
    Set();
  };

  using Sets = std::vector<Set>;

  /** the format is [MID[:INSERT:]]MSPEC
    where
    MID is a motif ID (optional)
    INSERT is list of comma separated numbers that give positions after which insert states are allowed
    MSPEC is one of three possibilities (tested in turn)
   * a path to a file. If such a file exists, it will be opened and a PWM will be read from it
   * a IUPAC regular expression
   * a length specification for lib Plasma, optionally specifying multiplicity
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
    enum class Kind {
      File,
      Seed,
      Plasma
    };
    Kind kind;
    std::string specification;
    std::string name;
    std::vector<size_t> insertions;
    std::vector<size_t> lengths;
    size_t multiplicity;
    Motif(const std::string &s="");
    Motif(const Motif &spec);
  };

  using Motifs = std::vector<Motif>;

  namespace Contrast {
    /**
      CONTRASTEXP is a contrast expression, as below
      CONTRASTEXP: ATOM | CONTRASTEXP-ATOM | CONTRASTEXP+ATOM
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

  /** the format is [MID[:CONTRASTEXP:]]OBJ
    where
    MID is a motif ID (optional)
    CONTRASTEXP is a contrast expression, as below
    OBJ is a string parsable as the desired objection function type (of template type X), i.e. it requires a definition of void parse_measure(const std::string&, &X)
    TODO: document that contrast names may not contain - or +
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

    Objective() :
      motif_name(""),
      contrast_expression(),
      measure(measure_t::Undefined)
    { };

    Objective(const std::string &s) :
      motif_name(""),
      contrast_expression(),
      measure()
    {
      std::string specification = s;
      size_t pos;
      if((pos = specification.find(":")) != std::string::npos) {
        motif_name = specification.substr(0,pos);
        specification = specification.substr(pos+1);
      }
      if((pos = specification.find(":")) != std::string::npos) {
        std::string expr_s = specification.substr(0,pos);
        contrast_expression = Contrast::make_contrast(expr_s);
        specification = specification.substr(pos+1);
      }
      parse_measure(specification, measure);
    };

    Objective(const Objective &obj) :
      motif_name(obj.motif_name),
      contrast_expression(obj.contrast_expression),
      measure(obj.measure)
    { };
  };

  template <typename X> typename Objective<X>::iterator begin(Objective<X> &objective) { return(begin(objective.contrast_expression)); }
  template <typename X> typename Objective<X>::iterator end(Objective<X> &objective) { return(end(objective.contrast_expression)); }
  template <typename X> typename Objective<X>::const_iterator begin(const Objective<X> &objective) { return(begin(objective.contrast_expression)); }
  template <typename X> typename Objective<X>::const_iterator end(const Objective<X> &objective) { return(end(objective.contrast_expression)); }

  std::istream &operator>>(std::istream &in, Motif &m);
  std::istream &operator>>(std::istream &in, Set &spec);
  template <typename X> std::istream &operator>>(std::istream &in, Objective<X> &spec) {
    std::string token;
    in >> token;
    spec = Objective<X>(token);
    return(in);
  }

  std::string to_string(const Motif &spec);
  std::string to_string(const Set &spec);
  template <typename X> std::string to_string(const Objective<X> &spec) {
    std::string s = std::string(spec.motif_name) + ":" + to_string(spec.contrast_expression) + ":" + measure2string(spec.measure);
    return(s);
  }

  std::ostream &operator<<(std::ostream &out, const Motif &m);
  std::ostream &operator<<(std::ostream &out, const Set &spec);
  template <typename X> std::ostream &operator<<(std::ostream &out, Objective<X> &spec) {
    out << to_string(spec);
    return(out);
  }

  /** A routine to interpret what the user has specified in terms of motifs, data sets, and objectives.
   * Contrasts:
   *   - find all contrast names occurring in data set specifications
   *   - give names to unnamed contrasts, making sure not to overwrite given ones
   *   - find all contrast names occurring in objective specifications
   *   - find all contrast names given in either data or objective specification
   *   - add all contrasts to those objectives that don't have any contrast annotated
   *   - check that the contrasts for the 2x2 measures are binary and that all discriminative ones are at least binary
   * Motifs:
   *   - get all motif names and check that none are duplicated
   *   - give names to unnamed motifs, making sure not to overwrite given ones
   *   - more that needs to be documented
   */
  template <typename X> void harmonize(Motifs &motifs, Sets &sets, std::vector<Objective<X>> &objectives, bool demultiplex, bool add_shuffles=true) {
    const bool debug = false;

    // find all contrast names occurring in data set specifications
    std::set<std::string> contrast_names_in_data;
    for(auto &spec: sets)
      if(spec.contrast != "")
        contrast_names_in_data.insert(spec.contrast);

    // give a name to unnamed contrasts, making sure not to overwrite given ones
    size_t contrast_idx = 0;
    std::string unnamed_name = "contrast" + boost::lexical_cast<std::string>(contrast_idx++);
    while(contrast_names_in_data.find(unnamed_name) != end(contrast_names_in_data))
      unnamed_name = "contrast" + boost::lexical_cast<std::string>(contrast_idx++);
    contrast_names_in_data.insert(unnamed_name);
    for(auto &spec: sets)
      if(spec.contrast == "")
        spec.contrast = unnamed_name;



    // find all contrast names occurring in objective specifications
    std::set<std::string> contrast_names_in_objectives;
    for(auto &objective: objectives)
      for(auto &atom: objective)
        contrast_names_in_objectives.insert(atom.contrast);

    if(debug) {
      std::cout << "Found contrast name in data specifications:";
      for(auto &s: contrast_names_in_data)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
      std::cout << "Found contrast name in objectives:";
      for(auto &s: contrast_names_in_objectives)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
    }


    // find all contrast names
    std::set<std::string> all_contrast_names;
    for(auto &s: contrast_names_in_data)
      all_contrast_names.insert(s);
    for(auto &s: contrast_names_in_objectives)
      all_contrast_names.insert(s);

    // add all contrasts to those objectives that don't have any contrast annotated
    for(auto &objective: objectives)
      if(objective.contrast_expression.empty())
        for(auto &s: all_contrast_names)
          objective.contrast_expression.push_back({+1, s});



    // check that the contrasts for the 2x2 measures are binary and that all discriminative ones are at least binary
    for(auto &objective: objectives)
      for(auto &atom: objective) {
        size_t contrast_size = 0;
        for(auto &spec: sets)
          if(atom.contrast == spec.contrast)
            contrast_size++;
        if(Measures::is_discriminative(objective.measure))
          if(contrast_size < 2) {
            if(add_shuffles and contrast_size == 1) {
              Set shuffle_spec;
              for(auto &spec: sets)
                if(atom.contrast == spec.contrast) {
                  shuffle_spec = Set(spec.path, true);
                  shuffle_spec.contrast = spec.contrast;
                  break;
                }
              sets.push_back(shuffle_spec);
            } else {
              std::cout << "Error: discriminative measure '" << objective.measure << "' requires at least a binary contrast in contrast '" << atom.contrast << "'." << std::endl;
              exit(-1);
            }
          }
        if(Measures::is_two_by_two(objective.measure))
          if(contrast_size != 2) {
            std::cout << "Error: 2x2 measure '" << objective.measure << "' only works on binary contrasts, and contrast '" << atom.contrast << "' is not binary." << std::endl;
            exit(-1);
          }
      }






    // get all motif names and check that none are duplicated
    std::set<std::string> motif_names_in_motifs;
    for(auto &motif: motifs)
      if(motif.name != "") {
        auto pair = motif_names_in_motifs.insert(motif.name);
        if(not pair.second) { // there was already a motif of that name
          std::cout << "Error: motif name not unique: '" << motif.name << "'." << std::endl;
          exit(-1);
        }
      }



    // give names to unnamed motifs, making sure not to overwrite given ones
    size_t motif_idx = 0;
    for(auto &motif: motifs)
      if(motif.name == "") {
        std::string name = "motif" + boost::lexical_cast<std::string>(motif_idx++);
        while(motif_names_in_motifs.find(name) != end(motif_names_in_motifs))
          name = "motif" + boost::lexical_cast<std::string>(motif_idx++);
        motif.name = name;
        motif_names_in_motifs.insert(name);
      }


    // check that no motifs have multiplicity zero
    for(auto &m: motifs)
      if(m.multiplicity == 0) {
        std::cout << "Error: multiplicity of motif \"" << m.name << "\" is zero!" << std::endl;
        exit(-1);
      }


    // if desired (=yes for Discrover, =no for Plasma) create multiple motifs for which multiplicity > 1
    std::map<std::string, std::list<std::string>> demux_motif_map;
    if(demultiplex) {
      Motifs new_motifs;
      // std::vector<Specification::Motif> new_motifs;
      for(auto &motif: motifs) {
        if(motif.multiplicity == 1)
          new_motifs.push_back(motif);
        else {
          for(size_t i = 0; i < motif.multiplicity; i++) {
            auto m = motif;
            m.multiplicity = 1;
            m.name += "_" + boost::lexical_cast<std::string>(i);

            new_motifs.push_back(m);

            demux_motif_map[motif.name].push_back(m.name);

            motif_names_in_motifs.insert(m.name);
            std::cout << "Demuxing: " << motif.name << " -> " << m.name << std::endl;
          }
          motif_names_in_motifs.erase(motif.name);
        }
      }
      motifs = new_motifs;
    }

    // TODO make sure that where motif names are mentioned they should be updated for the demultiplexed motif names

//
//      rg. objectives with empty contrast... they should be taken to mean all contrasts

    if(debug) {
      std::cout << "Found motif name in motif specification:";
      for(auto &s: motif_names_in_motifs)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
    }

//    if only -m mi is given then mi should be used for all motifs; i.e. objectives need to be constructed for each motif
//      otherwise, if there is -m "":mi and any other -m XYZ:mi statement it is an error
    if(objectives.size() == 1 and objectives[0].motif_name == "") {
      auto original_objective = *begin(objectives);
      objectives.clear();
      for(auto &name: motif_names_in_motifs) {
        Objective<X> objective(original_objective);
        objective.motif_name = name;
        objectives.push_back(objective);
      }
    }

    // find sequence sets that mention demuxed motif, and replace by new ones
    for(auto &spec: sets) {
      std::set<std::string> present_motifs;
      for(auto &m: spec.motifs)
        if(demux_motif_map.find(m) == end(demux_motif_map))
          present_motifs.insert(m);
        else
          for(auto &n: demux_motif_map.find(m)->second) {
            std::cout << "Demuxing (sets): " << m << " -> " << n << std::endl;
            present_motifs.insert(n);
          }
      spec.motifs = present_motifs;
    }

    // find objectives that mention demuxed motif, and replace by new ones
    if(demultiplex) {
      std::vector<Objective<X>> objs;
      for(auto &objective: objectives) {
        std::string motif_name = objective.motif_name;
        std::cout << "Check for Demuxing (obj): " << motif_name << std::endl;
        auto demux_motif_iter = demux_motif_map.find(motif_name);
        if(demux_motif_iter == end(demux_motif_map)) {
          objs.push_back(objective);
          std::cout << "Demuxing (obj): " << motif_name << " ok!" << std::endl;
        } else
          for(auto &demuxed: demux_motif_iter->second) {
            auto obj = objective;
            obj.motif_name = demuxed;
            std::cout << "Demuxing (obj): " << motif_name << " -> " << demuxed << std::endl;
            objs.push_back(obj);
          }
      }
      objectives = objs;
    }

    // check that no more than one objective exists for each motif, and that each objective refers to an existing motif
    // i.e. check that the relation between objective and motifs is one-to-one
    std::set<std::string> motif_names_in_objectives;
    for(auto &objective: objectives) {
      auto pair = motif_names_in_objectives.insert(objective.motif_name);
      if(not pair.second) { // there was already a motif of that name
        std::cout << "Error: motif name in objective not unique: '" << objective.motif_name << "'." << std::endl;
        exit(-1);
      }
      if(motif_names_in_motifs.find(objective.motif_name) == end(motif_names_in_motifs)) {
        std::cout << "Error: found a motif name in the objective for which no motif specification is found: '" << objective.motif_name << "'." << std::endl;
        exit(-1);
      }
    }
    if(motif_names_in_objectives.find("") != end(motif_names_in_objectives) and motif_names_in_objectives.size() > 1) {
      std::cout << "Error: when any objectives name motifs then all objectives must define motifs." << std::endl;
      exit(-1);
    }

    if(debug) {
      std::cout << "Found motif name in objectives:";
      for(auto &s: motif_names_in_objectives)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
    }
  }

};
#endif   /* ----- #ifndef SPECIFICATION_HPP ----- */

