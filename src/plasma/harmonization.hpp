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
 *       Filename:  harmonization.hpp
 *
 *    Description:  Ensure that motif, data, and objective specifications are consistent
 *
 *        Created:  Sun Jan 25 2015 16:05:34 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef HARMONIZATION_HPP
#define HARMONIZATION_HPP

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "specification.hpp"

namespace Specification {
namespace Exception {
namespace Discriminative {
template <typename measure_t>
struct UnaryContrast : public std::runtime_error {
  UnaryContrast(measure_t measure, const std::string &contrast)
      : runtime_error("Error: discriminative measure '"
                      + measure2string(measure)
                      + "' requires at least a binary contrast in contrast '"
                      + contrast + "'."){};
};
template <typename measure_t>
struct TwoByTwo : public std::runtime_error {
  TwoByTwo(measure_t measure, const std::string &contrast)
      : runtime_error("Error: 2x2 measure '" + measure2string(measure)
                      + "' only works on binary contrasts, and contrast '"
                      + contrast + "' is not binary."){};
};
}

namespace Motif {
struct NameNotUnique : public std::runtime_error {
  NameNotUnique(const std::string &name);
};
struct MultiplicityZero : public std::runtime_error {
  MultiplicityZero(const std::string &name);
};
struct NameNotUniqueInObjective : public std::runtime_error {
  NameNotUniqueInObjective(const std::string &name);
};
struct NoSpecfication : public std::runtime_error {
  NoSpecfication(const std::string &name);
};
}
namespace Objective {
struct MultpleObjectivesWithoutNamedMotifs : public std::runtime_error {
  MultpleObjectivesWithoutNamedMotifs();
};
struct NoMotif : public std::runtime_error {
  NoMotif(const std::string &token);
};
}
}

/** A routine to interpret what the user has specified in terms of motifs, data
 * sets, and objectives.
 * Contrasts:
 *   - find all contrast names occurring in data set specifications
 *   - create name for the unnamed contrast, ensuring not to overwrite given ones
 *   - find all contrast names occurring in objective specifications
 *   - find all contrast names given in either data or objective specification
 *   - add all contrasts to those objectives that don't have any contrast
 *     annotated
 *   - check that the contrasts for the 2x2 measures are binary and those for
 *     discriminative ones are at least binary
 * Motifs:
 *   - get all motif names and check that none are duplicated
 *   - give names to unnamed motifs, making sure not to overwrite given ones
 *   - check that no motifs have multiplicity zero
 *   - determine number of objectives that don't name motifs
 *     * if more than one such objectives exists, it is an error
 *     * if one such objective exists, then use it for all motifs for which no
 *       objective is given
 *   - check that no more than one objective exists for each motif, and that
 *     each objective refers to an existing motif
 *     i.e. check that the relation between objective and motifs is one-to-one
 *   - check at the end that no unnamed motifs are referenced in the objectives
 */
template <typename measure_t>
void harmonize(Motifs &motifs, Sets &sets,
               std::vector<Objective<measure_t>> &objectives,
               bool add_shuffles = true) {
  using namespace std;
  using Objectives  = vector<Objective<measure_t>>;
  const bool debug = false;

  // find all contrast names occurring in data set specifications
  set<string> contrast_names_in_data;
  bool found_unnamed_contrast_in_data = false;
  for (auto &spec : sets)
    if (spec.contrast != "")
      contrast_names_in_data.insert(spec.contrast);
    else
      found_unnamed_contrast_in_data = true;

  // create name for the unnamed contrast, ensuring not to overwrite given ones
  if (found_unnamed_contrast_in_data) {
    size_t contrast_idx = 0;
    string unnamed_name = "contrast" + std::to_string(contrast_idx++);
    while (contrast_names_in_data.find(unnamed_name)
           != end(contrast_names_in_data))
      unnamed_name = "contrast" + std::to_string(contrast_idx++);
    contrast_names_in_data.insert(unnamed_name);
    for (auto &spec : sets)
      if (spec.contrast == "")
        spec.contrast = unnamed_name;
  }

  // find all contrast names occurring in objective specifications
  set<string> contrast_names_in_objectives;
  for (auto &objective : objectives)
    for (auto &atom : objective)
      contrast_names_in_objectives.insert(atom.contrast);

  if (debug) {
    cout << "Found contrast name in data specifications:";
    for (auto &s : contrast_names_in_data)
      cout << " '" << s << "'";
    cout << endl;
    cout << "Found contrast name in objectives:";
    for (auto &s : contrast_names_in_objectives)
      cout << " '" << s << "'";
    cout << endl;
  }

  // find all contrast names given in either data or objective specification
  set<string> all_contrast_names;
  for (auto &s : contrast_names_in_data)
    all_contrast_names.insert(s);
  for (auto &s : contrast_names_in_objectives)
    all_contrast_names.insert(s);

  // add all contrasts to objectives that don't have any contrast annotated
  for (auto &objective : objectives)
    if (objective.contrast_expression.empty())
      for (auto &s : all_contrast_names)
        objective.contrast_expression.push_back({+1, s});

  // check:
  //   * contrasts for the 2x2 measures are binary
  //   * contrasts for discriminative measures are at least binary
  for (auto &objective : objectives)
    for (auto &atom : objective) {
      size_t contrast_size = 0;
      for (auto &spec : sets)
        if (atom.contrast == spec.contrast)
          contrast_size++;
      if (Measures::is_discriminative(objective.measure))
        if (contrast_size < 2) {
          if (add_shuffles and contrast_size == 1) {
            Set shuffle_spec;
            for (auto &spec : sets)
              if (atom.contrast == spec.contrast) {
                shuffle_spec = Set(spec.path, true);
                shuffle_spec.contrast = spec.contrast;
                break;
              }
            sets.push_back(shuffle_spec);
          } else
            throw Exception::Discriminative::UnaryContrast<measure_t>(
                objective.measure, atom.contrast);
        }
      if (Measures::is_two_by_two(objective.measure))
        if (contrast_size != 2)
          throw Exception::Discriminative::TwoByTwo<measure_t>(
              objective.measure, atom.contrast);
    }

  // get all motif names and check that none are duplicated
  set<string> motif_names_in_motifs;
  for (auto &motif : motifs)
    if (motif.name != "") {
      auto pair = motif_names_in_motifs.insert(motif.name);
      if (not pair.second)  // there was already a motif of that name
        throw Exception::Motif::NameNotUnique(motif.name);
    }

  // give names to unnamed motifs, making sure not to overwrite given ones
  size_t motif_idx = 0;
  for (auto &motif : motifs)
    if (motif.name == "") {
      string name = "motif" + std::to_string(motif_idx++);
      while (motif_names_in_motifs.find(name) != end(motif_names_in_motifs))
        name = "motif" + std::to_string(motif_idx++);
      motif.name = name;
      motif_names_in_motifs.insert(name);
    }

  // check that no motifs have multiplicity zero
  for (auto &m : motifs)
    if (m.multiplicity == 0)
      throw Exception::Motif::MultiplicityZero(m.name);

  // TODO? rg. objectives with empty contrast: they are taken to mean all contrasts

  if (debug) {
    cout << "Found motif name in motif specification:";
    for (auto &s : motif_names_in_motifs)
      cout << " '" << s << "'";
    cout << endl;
  }

  // determine number of objectives that don't name motifs
  // If an objective exists that does not name motifs, then use this for all
  // motifs for which no objective is given
  Objectives objectives_without_named_motifs;
  vector<string> motifs_in_objectives;
  for (auto &objective : objectives)
    if (objective.motif_name == "")
      objectives_without_named_motifs.push_back(objective);
    else
      motifs_in_objectives.push_back(objective.motif_name);

  switch (objectives_without_named_motifs.size()) {
    case 0:
      break;
    case 1:
      // use this objective for all motifs for which no objective is given
      {
        auto obj_wo_named_motfs = *begin(objectives_without_named_motifs);
        // copy objectives
        auto original_objectives = objectives;
        // reset objectives
        objectives.clear();
        // add objectives with named motifs
        for (auto &objective : objectives)
          if (objective.motif_name != "")
            objectives.push_back(objective);
        // generate an objective for each motif not named in an objective
        for (auto &name : motif_names_in_motifs)
          if (find(begin(motifs_in_objectives), end(motifs_in_objectives), name)
              == end(motifs_in_objectives)) {
            auto objective = obj_wo_named_motfs;
            objective.motif_name = name;
            objectives.push_back(objective);
          }
      }
      break;
    default:
      throw Exception::Objective::MultpleObjectivesWithoutNamedMotifs();
  }

  // check that no more than one objective exists for each motif, and that each
  // objective refers to an existing motif
  // i.e. check that the relation between objective and motifs is one-to-one
  set<string> motif_names_in_objectives;
  for (auto &objective : objectives) {
    auto pair = motif_names_in_objectives.insert(objective.motif_name);
    if (not pair.second)  // there was already a motif of that name
      throw Exception::Motif::NameNotUniqueInObjective(objective.motif_name);
    if (motif_names_in_motifs.find(objective.motif_name)
        == end(motif_names_in_motifs))
      throw Exception::Motif::NoSpecfication(objective.motif_name);
  }

  // finally check that no unnamed motifs are referenced in the objectives
  for (auto &objective : objectives)
    if (objective.motif_name == "")
      throw Exception::Objective::NoMotif(to_string(objective));

  if (debug) {
    cout << "Found motif name in objectives:";
    for (auto &s : motif_names_in_objectives)
      cout << " '" << s << "'";
    cout << endl;
  }
}
}

#endif /* ----- #ifndef HARMONIZATION_HPP ----- */
