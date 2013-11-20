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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef SPECIFICATION_HPP
#define SPECIFICATION_HPP

#include <iostream> // TODO remove
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <set>
#include <map>

namespace Measures {
  template <typename X> bool is_discriminative(X measure);
  template <typename X> bool is_two_by_two(X measure);
};

namespace Specification {
  /** the format is [MIDs[:SERIES:]]PATH
   * where
   * MIDs a set of motif ID (optional)
   * SERIES is the name of a series this data set belongs to; if unspecified, the unnamed general series is used
   * PATH a path to a FASTA file
   */
  struct DataSet {
    std::string series;
    std::string path;
    bool is_shuffle;
    std::set<std::string> motifs;
    DataSet(const std::string &s, bool shuffled=false);
    DataSet(const DataSet &spec);
    DataSet();
  };

  typedef std::vector<DataSet> DataSets;

  /** the format is [MID[:INSERT:]]MSPEC
    where
    MID is a motif ID (optional)
    INSERT is list of comma separated numbers that give positions after which insert states are allowed
    MSPEC is one of three possibilities (tested in turn)
   * a path to a file. If such a file exists, it will be opened and a PWM will be read from it
   * a IUPAC regular expression
   * a length specification for lib Plasma
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
    Motif(const std::string &s="");
    Motif(const Motif &spec);
  };

  typedef std::vector<Motif> Motifs;

  namespace Series {
    /**
      SERIESEXP is a series expression, as below
      SERIESEXP: ATOM | SERIESEXP-ATOM | SERIESEXP+ATOM
     */
    struct Atom {
      int sign;
      std::string series;
    };

    typedef std::vector<Atom> Expression;
    Expression make_series(const std::string &path);

    std::string to_string(const Atom &atom);
    std::string to_string(const Expression &expr);
    std::ostream &operator<<(std::ostream &out, const Atom &atom);
    std::ostream &operator<<(std::ostream &out, const Expression &expr);
  }

  /** the format is [MID[:SERIESEXP:]]OBJ
    where
    MID is a motif ID (optional)
    SERIESEXP is a series expression, as below
    OBJ is a string parsable as the desired objection function type (of template type X), i.e. it requires a definition of void parse_measure(const std::string&, &X)
    TODO: document that series names may not contain - or +
   */
  template <typename X>
  struct Objective {
    typedef X measure_t;
    typedef Series::Atom series_atom_t;
    typedef Series::Expression series_expression_t;
    typedef typename series_expression_t::iterator iterator;
    typedef typename series_expression_t::const_iterator const_iterator;

    std::string motif_name;
    series_expression_t series_expression;
    measure_t measure;

    Objective() :
      motif_name(""),
      series_expression(),
      measure(measure_t::Undefined)
    { };

    Objective(const std::string &s) :
      motif_name(""),
      series_expression(),
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
        series_expression = Series::make_series(expr_s);
        specification = specification.substr(pos+1);
      }
      parse_measure(specification, measure);
    };

    Objective(const Objective &obj) :
      motif_name(obj.motif_name),
      series_expression(obj.series_expression),
      measure(obj.measure)
    { };
  };

  template <typename X> typename Objective<X>::iterator begin(Objective<X> &objective) { return(begin(objective.series_expression)); }
  template <typename X> typename Objective<X>::iterator end(Objective<X> &objective) { return(end(objective.series_expression)); }
  template <typename X> typename Objective<X>::const_iterator begin(const Objective<X> &objective) { return(begin(objective.series_expression)); }
  template <typename X> typename Objective<X>::const_iterator end(const Objective<X> &objective) { return(end(objective.series_expression)); }

  std::istream &operator>>(std::istream &in, Motif &m);
  std::istream &operator>>(std::istream &in, DataSet &spec);
  template <typename X> std::istream &operator>>(std::istream &in, Objective<X> &spec) {
    std::string token;
    in >> token;
    spec = Objective<X>(token);
    return(in);
  }

  std::string to_string(const Motif &spec);
  std::string to_string(const DataSet &spec);
  template <typename X> std::string to_string(const Objective<X> &spec) {
    std::string s = std::string(spec.motif_name) + ":" + to_string(spec.series_expression) + ":" + measure2string(spec.measure);
    return(s);
  }

  std::ostream &operator<<(std::ostream &out, const Motif &m);
  std::ostream &operator<<(std::ostream &out, const DataSet &spec);
  template <typename X> std::ostream &operator<<(std::ostream &out, Objective<X> &spec) {
    out << to_string(spec);
    return(out);
  }

  /** A routine to interpret what the user has specified in terms of motifs, data sets, and objectives.
   * Series:
   *   - find all series names occurring in data set specifications
   *   - give names to unnamed series, making sure not to overwrite given ones
   *   - find all series names occurring in objective specifications
   *   - find all series names given in either data or objective specification
   *   - add all series to those objectives that don't have any series annotated
   *   - check that the series for the 2x2 measures are binary and that all discriminative ones are at least binary
   * Motifs:
   *   - get all motif names and check that none are duplicated
   *   - give names to unnamed motifs, making sure not to overwrite given ones
   *   - more that needs to be documented
   */
  template <typename X> void harmonize(Motifs &motifs, DataSets &data, std::vector<Objective<X>> &objectives, bool add_shuffles=true) {
    const bool debug = false;

    // find all series names occurring in data set specifications
    std::set<std::string> series_names_in_data;
    for(auto &spec: data)
      if(spec.series != "")
        series_names_in_data.insert(spec.series);

    // give a name to unnamed series, making sure not to overwrite given ones
    size_t contrast_idx = 0;
    std::string unnamed_name = "contrast" + boost::lexical_cast<std::string>(contrast_idx++);
    while(series_names_in_data.find(unnamed_name) != end(series_names_in_data))
      unnamed_name = "contrast" + boost::lexical_cast<std::string>(contrast_idx++);
    series_names_in_data.insert(unnamed_name);
    for(auto &spec: data)
      if(spec.series == "")
        spec.series = unnamed_name;



    // find all series names occurring in objective specifications
    std::set<std::string> series_names_in_objectives;
    for(auto &objective: objectives)
      for(auto &atom: objective)
        series_names_in_objectives.insert(atom.series);

    if(debug) {
      std::cout << "Found series name in data specifications:";
      for(auto &s: series_names_in_data)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
      std::cout << "Found series name in objectives:";
      for(auto &s: series_names_in_objectives)
        std::cout << " '" << s << "'";
      std::cout << std::endl;
    }


    // find all series names
    std::set<std::string> all_series_names;
    for(auto &s: series_names_in_data)
      all_series_names.insert(s);
    for(auto &s: series_names_in_objectives)
      all_series_names.insert(s);

    // add all series to those objectives that don't have any series annotated
    for(auto &objective: objectives)
      if(objective.series_expression.empty())
        for(auto &s: all_series_names)
          objective.series_expression.push_back({+1, s});



    // check that the series for the 2x2 measures are binary and that all discriminative ones are at least binary
    for(auto &objective: objectives)
      for(auto &atom: objective) {
        size_t series_size = 0;
        for(auto &spec: data)
          if(atom.series == spec.series)
            series_size++;
        if(Measures::is_discriminative(objective.measure))
          if(series_size < 2) {
            if(add_shuffles and series_size == 1) {
              DataSet shuffle_spec;
              for(auto &spec: data)
                if(atom.series == spec.series) {
                  shuffle_spec = DataSet(spec.path, true);
                  shuffle_spec.series = spec.series;
                  break;
                }
              data.push_back(shuffle_spec);
            } else {
              std::cout << "Error: discriminative measure '" << objective.measure << "' requires at least a binary contrast in series '" << atom.series << "'." << std::endl;
              exit(-1);
            }
          }
        if(Measures::is_two_by_two(objective.measure))
          if(series_size != 2) {
            std::cout << "Error: 2x2 measure '" << objective.measure << "' only works on binary contrasts, and series '" << atom.series << "' is not binary." << std::endl;
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
    for(auto &motif:motifs)
      if(motif.name == "") {
        std::string name = "motif" + boost::lexical_cast<std::string>(motif_idx++);
        while(motif_names_in_motifs.find(name) != end(motif_names_in_motifs))
          name = "motif" + boost::lexical_cast<std::string>(motif_idx++);
        motif.name = name;
        motif_names_in_motifs.insert(name);
      }




//
//      rg. objectives with empty series... they should be taken to mean all series

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

