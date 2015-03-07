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
 *       Filename:  hmm_aux.cpp
 *
 *    Description:  Auxiliary routines of the HMM class
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iomanip>
#include <set>
#include "../aux.hpp"
#include "hmm.hpp"

using namespace std;

bitmask_t HMM::compute_bitmask(const Training::Task &task) const {
  bitmask_t present = 0;
  for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if (task.motif_name == groups[group_idx].name)
      present[group_idx] = 1;
  return present;
}

void HMM::shift_forward(size_t group_idx, size_t n) {
  for (size_t i = 0; i < groups[group_idx].states.size() - 1; i++)
    for (size_t j = 0; j < n_emissions; j++)
      emission(groups[group_idx].states[i], j)
          = emission(groups[group_idx].states[i + 1], j);
  for (size_t j = 0; j < n_emissions; j++)
    emission(*groups[group_idx].states.rbegin(), j) = 1.0 / n_emissions;
  normalize_emission(emission);
}

void HMM::shift_backward(size_t group_idx, size_t n) {
  for (size_t i = groups[group_idx].states.size() - 1; i > 0; i--)
    for (size_t j = 0; j < n_emissions; j++)
      emission(groups[group_idx].states[i - 1], j)
          = emission(groups[group_idx].states[i], j);
  for (size_t j = 0; j < n_emissions; j++)
    emission(*groups[group_idx].states.begin(), j) = 1.0 / n_emissions;
  normalize_emission(emission);
}

void HMM::serialize(ostream &os, const ExecutionInformation &exec_info,
                    size_t format_version) const {
  const ios::fmtflags flags(os.flags());
  const string param_format_string = "# HMM parameter format version ";
  const string transition_matrix_string = "Transition matrix";
  const string emission_matrix_string = "Emission matrix";
  const streamsize output_precision = 15;
  switch (format_version) {
    // TODO bump version because no order is stored?
    case 6:
      os << param_format_string << format_version << endl;
      os << "# " << exec_info.program_name << " " << exec_info.hmm_version
         << " [" << exec_info.git_branch << " branch]" << endl;
      os << "# Run on " << exec_info.datetime << endl;
      os << "# Run in " << exec_info.directory << endl;
      os << "# Command = " << exec_info.cmdline << endl;
      for (auto &x : registration.datasets) {
        os << "Dataset " << x.second.spec.path << " " << x.first
           << " class = " << x.second.class_prior << " motif = ";
        for (auto &y : x.second.motif_prior)
          if (y.first != 0)
            os << " " << y.first << "/" << y.second;
        os << endl;
      }
      os << n_states << " states" << endl;
      os << n_emissions << " emissions" << endl;
      for (auto &motif : groups) {
        os << "Motif \"" << motif.name << "\"";
        for (auto &state : motif.states)
          os << " " << state;
        os << endl;
      }
      os << "State class";
      for (auto &motif : group_ids)
        os << " " << motif;
      os << endl;
      os << transition_matrix_string << endl;
      os.precision(output_precision);
      for (size_t i = 0; i < n_states; i++) {
        os << setw(3) << fixed << i;
        for (size_t j = 0; j < n_states; j++)
          os << " " << setw(2 + output_precision) << fixed << transition(i, j);
        os << endl;
      }
      os << emission_matrix_string << endl;
      for (size_t i = 0; i < n_states; i++) {
        os << setw(3) << fixed << i;
        for (size_t j = 0; j < n_emissions; j++)
          os << " " << setw(2 + output_precision) << fixed << emission(i, j);
        os << endl;
      }
      break;
    default:
      throw Exception::HMM::ParameterFile::UnsupportedVersion(format_version);
  }
  os.flags(flags);
};

void HMM::deserialize(istream &is) {
  groups = vector<Group>();
  group_ids = vector<size_t>();

  const string param_format_string = "# HMM parameter format version ";
  string line;
  safeGetline(is, line);
  if (line.length() > param_format_string.length()
      and line.substr(0, param_format_string.length()) == param_format_string) {
    size_t format_version
        = atoi(line.substr(param_format_string.length()).c_str());
    if (format_version < 4)
      throw Exception::HMM::ParameterFile::UnsupportedVersion(format_version);
    while (is.good() and line[0] == '#')
      safeGetline(is, line);

    if (format_version >= 5) {
      const string key_a = "Class_prior ";
      const string key_b = "Motif_prior[0] ";
      const string key_c = "Motif_prior[1] ";
      // TODO parse information regarding registered data sets
      while (line.find(" states") == string::npos)
        safeGetline(is, line);
      if (false) {
        // if(line.substr(0, key_a.size()) == key_a)
        //   class_prior = atof(line.substr(key_a.size()).c_str());
        // else
        //   cout << "Error trying to parse HMM parameters: class prior could
        //   not be parsed: '" << line << "'." << endl;
        safeGetline(is, line);
        // TODO parse information regarding registered data sets
        // if(line.substr(0, key_b.size()) == key_b)
        //   motif_prior[0] = atof(line.substr(key_b.size()).c_str());
        // else
        //   cout << "Error trying to parse HMM parameters: first conditional
        //   motif prior could not be parsed: '" << line << "'." << endl;
        safeGetline(is, line);
        // TODO parse information regarding registered data sets
        // if(line.substr(0, key_c.size()) == key_c)
        //   motif_prior[1] = atof(line.substr(key_c.size()).c_str());
        // else
        //   cout << "Error trying to parse HMM parameters: second conditional
        //   motif prior could not be parsed: '" << line << "'." << endl;
        safeGetline(is, line);
      }
    }

    if (verbosity >= Verbosity::debug)
      cerr << "n_states line = " << line << endl;
    n_states = atoi(line.c_str());
    if (verbosity >= Verbosity::debug)
      cerr << "n_states = " << n_states << endl;
    safeGetline(is, line);
    if (verbosity >= Verbosity::debug)
      cerr << "line = " << line << endl;
    size_t n_emis = atoi(line.c_str());
    if (n_emis != n_emissions)
      throw Exception::HMM::ParameterFile::SyntaxError(
          "this version only works with " + to_string(n_emissions)
          + " emissions, while the .hmm file specifies " + to_string(n_emis)
          + ".");

    if (verbosity >= Verbosity::debug)
      cout << "n_states = " << n_states << " n_emissions " << n_emissions
           << endl;

    safeGetline(is, line);
    while (line.substr(0, 5) == "Motif") {
      size_t start = line.find('"');
      size_t stop = line.find('"', start + 1);
      if (verbosity >= Verbosity::debug)
        cout << "start = " << start << " stop " << stop << endl;
      string name = line.substr(start + 1, stop - start - 1);
      if (verbosity >= Verbosity::debug)
        cout << "name = '" << name << "' line " << line << endl;
      line = line.substr(stop + 2);
      if (verbosity >= Verbosity::debug)
        cout << "name = " << name << " line " << line << endl;
      vector<size_t> states;
      while (line != "") {
        start = line.find(' ');
        if (verbosity >= Verbosity::debug)
          cout << "start = " << start << " line " << line << endl;
        if (start > 0) {
          size_t state = atoi(line.substr(0, start).c_str());
          states.push_back(state);
        }
        if (start == string::npos)
          break;
        else
          line = line.substr(start + 1);
      }
      Group::Kind kind = Group::Kind::Motif;
      if (name == "Special")
        kind = Group::Kind::Special;
      else if (name == "Background")
        kind = Group::Kind::Background;
      Group m = {kind, name, states};
      groups.push_back(m);
      safeGetline(is, line);
    }

    // safeGetline(is, line);
    if (line.substr(0, 11) != "State class")
      throw Exception::HMM::ParameterFile::SyntaxError(
          "expected \"State class\" but found instead: \"" + line + "\".");

    line = line.substr(12);
    while (line != "") {
      size_t start = line.find(" ");
      if (start > 0) {
        size_t motif = atoi(line.substr(0, start).c_str());
        group_ids.push_back(motif);
      }
      if (start == string::npos)
        break;
      else
        line = line.substr(start + 1);
    }

    safeGetline(is, line);
    if (line.substr(0, 11) == "State order") {
      cout << "Found a state order line; checking that all orders are zero, "
              "and ignoring." << endl;

      line = line.substr(12);
      while (line != "") {
        size_t start = line.find(" ");
        if (start > 0) {
          size_t o = atoi(line.substr(0, start).c_str());
          if (o != 0)
            throw Exception::HMM::ParameterFile::SyntaxError(
                "found a non-zero state order. This version version only works "
                "for zeroth order.");
        }
        if (start == string::npos)
          break;
        else
          line = line.substr(start + 1);
      }
      safeGetline(is, line);
    }

    last_state = n_states - 1;

    transition.resize(n_states, n_states);
    emission.resize(n_states, n_emissions);

    if (line != "Transition matrix")
      throw Exception::HMM::ParameterFile::SyntaxError(
                                 "expecting line \"Transition matrix\", "
                                 "but got instead \"" + line + "\".");
    for (size_t i = 0; i < n_states; i++) {
      size_t tmp;
      is >> tmp;
      for (size_t j = 0; j < n_states; j++)
        is >> transition(i, j);
    }
    safeGetline(is, line);
    safeGetline(is, line);
    if (line != "Emission matrix")
      throw Exception::HMM::ParameterFile::SyntaxError(
          "expecting line \"Emission matrix\", but got instead \"" + line
          + "\".");
    for (size_t i = 0; i < n_states; i++) {
      size_t tmp;
      is >> tmp;
      for (size_t j = 0; j < n_emissions; j++)
        is >> emission(i, j);
    }
  } else {
    is >> transition;
    is >> emission;
    last_state = transition.size1() - 1;
  }
  if (n_states != transition.size1())
      throw Exception::HMM::ParameterFile::SyntaxError("inconsistent number of states.");

  finalize_initialization();
};

string HMM::path2string_state(const HMM::StatePath &path) const {
  const char start_symb = '^';
  const char bg_symb = '0';
  string p;
  for (size_t j = 0; j < path.size(); j++)
    switch (path[j]) {
      case start_state:
        p += start_symb;
        break;
      case bg_state:
        p += bg_symb;
        break;
      default:
        char x = path[j] - first_state;
        char c = x + 'A';
        if (x > 'Z') {
          c = x - ('Z' - 'A' + 1) + 'a';
          if (c > 'z') {
            cout << "Error: state index is too large" << endl;
            c = '*';
          }
        }
        p += c;
        break;
    };
  return p;
};

string HMM::path2string_group(const HMM::StatePath &path) const {
  const char start_symb = '^';
  const char bg_symb = '-';
  char first_symb = '0';
  if (groups.size() - 2 > 10)
    first_symb = 'A';
  string p;
  for (size_t j = 0; j < path.size(); j++)
    switch (path[j]) {
      case start_state:
        p += start_symb;
        break;
      case bg_state:
        p += bg_symb;
        break;
      default:
        char c = first_symb + (group_ids[path[j]] - 2);
        p += c;
        break;
    };
  return p;
};

void HMM::to_dot(ostream &os, double minimum_transition) const {
  static const string digits = "0123456789";
  os << "digraph d {" << endl;
  os << "  rankdir=\"LR\";" << endl;
  for (size_t i = 0; i < n_states; i++)
    for (size_t j = 0; j < n_states; j++) {
      if (transition(i, j) >= minimum_transition) {
        string label1;
        switch (i) {
          case start_state:
            label1 = "S";
            break;
          default:
            if (i < first_state)
              label1 = string("B") + digits[i - bg_state];
            else
              label1 = string() + digits[i - first_state];
            break;
        }
        string label2;
        switch (j) {
          case start_state:
            label2 = "S";
            break;
          default:
            if (j < first_state)
              label2 = string("B") + digits[j - bg_state];
            else
              label2 = string() + digits[j - first_state];
            break;
        }
        os << label1 << " -> " << label2 << ";" << endl;
      }
    }
  os << "}" << endl;
};

bool HMM::check_consistency_transitions(double eps) const {
  for (size_t k = 0; k < n_states; k++) {
    double p = 0;
    for (size_t l = 0; l < n_states; l++)
      p += transition(k, l);
    if (fabs(1 - p) > eps)
      return false;
  }
  return true;
}

bool HMM::check_consistency_emissions(double eps) const {
  for (size_t k = 0; k < n_states; k++) {
    double p = 0;
    for (size_t l = 0; l < n_emissions; l++)
      p += emission(k, l);
    if (fabs(1 - p) > eps)
      return false;
  }
  return true;
}

bool HMM::check_consistency(double eps) const {
  bool ok = true;
  ok = ok and check_consistency_transitions(eps);
  ok = ok and check_consistency_emissions(eps);
  ok = ok and group_ids.size() == n_states;
  for (auto &m : groups)
    for (auto &state : m.states)
      if (state > n_states)
        ok = false;
  for (auto &d : group_ids)
    if (d > groups.size())
      ok = false;
  return ok;
}

ostream &operator<<(ostream &os, const HMM &hmm) {
  os << "HMM: " << hmm.n_states << " states (including start state), "
     << hmm.n_emissions << " emissions." << endl
     << "Transition probability matrix: " << hmm.transition << endl
     << "Emission probability matrix: " << hmm.emission << endl;
  return os;
}

/** Set the emissions according to the desired matrix.
 * Also embed the matrix with pad_left uniform emission positions on the left,
 * and pad_left uniform emission positions on the right. Finally, a number of
 * n_insertions positions with uniform emissions are added for insertion
 * states.
 */
void HMM::set_motif_emissions(const matrix_t &e, size_t first_padded,
                              size_t n_insertions, size_t pad_left,
                              size_t pad_right) {
  if (verbosity >= Verbosity::debug)
    cout << "set motif_emissions; n_insertions = " << n_insertions << endl;

  const size_t n = e.size1();

  // first set the emissions of the padding states on the left to uniform
  for (size_t i = 0; i < pad_left; i++) {
    for (size_t j = 0; j < e.size2(); j++)
      emission(first_padded + i, j) = 1.0 / e.size2();
    for (size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + i, j) = 0;
  }

  // set the emissions of the states of the actual motif to the given matrix
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + i, j) = e(i, j);
    for (size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + i, j) = 0;
  }

  // finally set the emissions of the padding states on the right to uniform
  for (size_t i = 0; i < pad_right; i++) {
    for (size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + n + i, j) = 1.0 / e.size2();
    for (size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + n + i, j) = 0;
  }

  // set the emissions of the insert states to uniform
  for (size_t i = 0; i < n_insertions; i++) {
    for (size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + n + pad_right + i, j) = 1.0
                                                                 / e.size2();
    for (size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + n + pad_right + i, j) = 0;
  }
}

void HMM::normalize_transition(matrix_t &m) const {
  for (size_t k = 0; k < m.size1(); k++) {
    double z = 0;
    for (size_t l = 0; l < m.size2(); l++)
      z += m(k, l);
    if (z > 0)
      for (size_t l = 0; l < m.size2(); l++)
        m(k, l) = m(k, l) / z;
    else
      for (size_t l = 0; l < m.size2(); l++)
        m(k, l) = 0;
  }
}

void HMM::normalize_emission(matrix_t &m) const {
  for (size_t k = bg_state; k < m.size1(); k++) {
    if (verbosity >= Verbosity::debug)
      cout << "k = " << k << endl;
    double z = 0;
    for (size_t b = 0; b < n_emissions; b++)
      z += m(k, b);
    if (z > 0)
      for (size_t b = 0; b < n_emissions; b++)
        m(k, b) /= z;
    else
      for (size_t b = 0; b < n_emissions; b++)
        m(k, b) = 0;
  }

  // enforce emissions of start state
  for (size_t b = 0; b < n_emissions; b++)
    m(start_state, b) = 0;
}

Training::Tasks HMM::define_training_tasks(const Options::HMM &options) const {
  Training::Tasks tasks;

  // First the (mostly) discriminative tasks
  bool atleast_one_discriminative_task = false;
  bool do_transition = options.bg_learning == Training::Method::Gradient
                       and options.objectives.size() == 1;
  for (auto &objective : options.objectives) {
    if (not Measures::is_discriminative(objective.measure))
      continue;
    atleast_one_discriminative_task = true;
    Training::Task task;
    task.motif_name = objective.motif_name;
    task.measure = objective.measure;
    task.contrast_expression = objective.contrast_expression;
    bool found = false;
    for (auto &group : groups) {
      if (group.name == objective.motif_name) {
        found = true;
        copy(begin(group.states), end(group.states),
             back_inserter(task.targets.emission));
        if (do_transition)
          copy(begin(group.states), end(group.states),
               back_inserter(task.targets.transition));
      }
    }

    if (not found) {
      if (options.verbosity >= Verbosity::verbose)
        cout << "Skipping discriminative training task for motif "
             << task.motif_name << " because it is not represented in the HMM."
             << endl;
    } else {
      tasks.push_back(task);

      if (options.verbosity >= Verbosity::verbose) {
        cout << "Generated primary training targets." << endl << "Emissions:";
        for (auto e : task.targets.emission)
          cout << " " << e;
        cout << endl;
        cout << "Transitions:";
        for (auto t : task.targets.transition)
          cout << " " << t;
        cout << endl;
      }
    }
  }

  // get all assigned emission and transition states
  // check that none are doubly assigned
  set<size_t> e_states, t_states;
  for (auto &task : tasks) {
    for (auto &e_state : task.targets.emission) {
      auto pair = e_states.insert(e_state);
      if (not pair.second)
        throw Exception::HMM::Learning::MultipleTasks("emission");
    }
    for (auto &t_state : task.targets.transition) {
      auto pair = t_states.insert(t_state);
      if (not pair.second)
        throw Exception::HMM::Learning::MultipleTasks("transition");
    }
  }

  // get all contrast names
  set<string> contrast_names;
  for (auto &task : tasks)
    for (auto &expr : task)
      contrast_names.insert(expr.contrast);

  // And then the generative part
  if (options.bg_learning != Training::Method::None) {
    Training::Task task;
    if (atleast_one_discriminative_task or options.objectives.empty()) {
      task.motif_name = "Context";
      task.measure = Measure::Likelihood;
    } else {
      task.motif_name = options.objectives[0].motif_name;
      task.measure = options.objectives[0].measure;
    }

    for (auto &contrast_name : contrast_names)
      task.contrast_expression.push_back({+1, contrast_name});
    for (size_t i = 0; i < n_states; i++) {
      if (e_states.find(i) == end(e_states))
        if (i > start_state)
          task.targets.emission.push_back(i);
      if (t_states.find(i) == end(t_states))
        task.targets.transition.push_back(i);
    }

    tasks.push_back(task);

    if (options.verbosity >= Verbosity::verbose) {
      cout << "Generated secondary training targets." << endl << "Emissions:";
      for (auto e : task.targets.emission)
        cout << " " << e;
      cout << endl;
      cout << "Transitions:";
      for (auto t : task.targets.transition)
        cout << " " << t;
      cout << endl;
    }
  }
  return tasks;
}

size_t HMM::non_zero_parameters(
    const Training::Targets &training_targets) const {
  size_t z = 0;
  for (auto i : training_targets.transition)
    for (size_t j = 0; j < n_states; j++)
      if (transition(i, j) > 0)
        z++;
  for (auto i : training_targets.emission)
    for (size_t j = 0; j < n_emissions; j++)
      if (emission(i, j) > 0)
        z++;
  if (verbosity >= Verbosity::verbose)
    cout << "non_zero_parameters = " << z << endl;
  return z;
}

size_t HMM::n_parameters() const {
  size_t n = n_states * (n_states - 1) + n_states * (n_emissions - 1);
  if (verbosity >= Verbosity::verbose)
    cout << "n_parameters = " << n << endl;
  return n;
}

size_t HMM::count_motif(const HMM::StatePath &path, size_t motif) const {
  if (motif >= groups.size())
    return 0;
  size_t first_state
      = *groups[motif].states.begin();  // TODO make this more flexible
  size_t n = 0;
  for (auto x : path)
    if (x == first_state)
      n++;
  return n;
}

bool HMM::is_motif_group(size_t group_idx) const {
  return groups[group_idx].kind == Group::Kind::Motif;
}

bool HMM::is_motif_state(size_t state_idx) const {
  return is_motif_group(group_ids[state_idx]);
}

HMM::mask_t HMM::compute_mask(const Data::Collection &collection) const {
  if (verbosity >= Verbosity::debug)
    cout << "HMM::compute-mask(Data::Collection)" << endl;
  HMM::mask_t mask;
  for (auto &contrast : collection)
    for (auto &dataset : contrast) {
      mask_sub_t m;
      for (auto &seq : dataset) {
        vector<size_t> v;
        HMM::StatePath path;
        viterbi(seq, path);
        size_t idx = 0;
        for (auto state : path) {
          if (state >= first_state)
            v.push_back(idx);  // TODO: check: should I insert idx - 1 instead?
          // -> No; but there was an issue with reverse complements; fixed in
          // src/plasma/data.cpp DataSet::mask by pos = 2 * n - pos where n is
          // the sequence length
          idx++;
        }
        if (not v.empty())
          m[seq.definition] = v;
      }
      if (not m.empty())
        mask[dataset.path] = m;
    }
  return mask;
}

Training::Range HMM::complementary_states(size_t group_idx) const {
  Training::Range range;
  for (size_t i = start_state; i < n_states; i++)
    if (not is_motif_state(i) or group_ids[i] != group_idx)
      range.push_back(i);
  if (verbosity >= Verbosity::debug) {
    cout << "Complementary states =";
    for (auto &x : range)
      cout << " " << x;
    cout << endl;
  }
  return range;
}

Training::Range HMM::complementary_states_mask(bitmask_t present_mask) const {
  Training::Range range;
  for (size_t i = start_state; i < n_states; i++)
    range.push_back(i);
  auto in_mask = [&](size_t state) {
    return (present_mask & bitmask_t(1 << group_ids[state])) != 0;
  };
  range.erase(remove_if(begin(range), end(range), in_mask), end(range));
  if (verbosity >= Verbosity::debug) {
    cout << "Complementary states =";
    for (auto &x : range)
      cout << " " << x;
    cout << endl;
  }
  return range;
}

// Find those states reachable from the set of states given as argument.
// NOTE: There's a slight hack in that transitions from the set of final
// states to the set of initial states are not considered.
map<size_t, set<size_t>> reachable_states(const matrix_t &transition,
                                          const vector<size_t> &states,
                                          const set<size_t> &initial,
                                          const set<size_t> &final) {
  map<size_t, set<size_t>> r;
  for (auto &i : states)
    for (auto &j : states)
      if (transition(i, j) > 0)
        if (final.find(i) == end(final) or initial.find(j) == end(initial))
          r[i].insert(j);
  size_t n = states.size();
  while (n--) {
    for (auto &p : r) {
      set<size_t> y;
      for (auto &z : p.second) {
        y.insert(z);
        for (auto &a : r[z])
          y.insert(a);
      }
      p.second = y;
    }
  }
  return r;
}

set<size_t> reachable_from_state(const matrix_t &transition,
                                 const size_t &state) {
  set<size_t> r;
  for (size_t i = 0; i < transition.size2(); i++)
    if (transition(state, i) > 0)
      r.insert(i);
  return r;
}

set<size_t> initial_states(const matrix_t &transition,
                           const vector<size_t> &states) {
  const size_t bg_state = 1;  // HMM::bg_state is protected
  auto reachable_from_bg = reachable_from_state(transition, bg_state);
  set<size_t> initial;
  for (auto state : states)
    if (reachable_from_bg.find(state) != end(reachable_from_bg))
      initial.insert(state);
  return initial;
}

set<size_t> final_states(const matrix_t &transition,
                         const vector<size_t> &states) {
  const size_t bg_state = 1;  // HMM::bg_state is protected
  set<size_t> final;
  for (auto state : states)
    if (transition(state, bg_state) > 0)
      final.insert(state);
  return final;
}

vector<size_t> topological_order(const matrix_t &transition,
                                 const vector<size_t> &states) {
  auto initial = initial_states(transition, states);
  auto final = final_states(transition, states);
  auto reachable = reachable_states(transition, states, initial, final);
  vector<size_t> order = states;
  sort(begin(order), end(order), [&reachable](size_t x, size_t y) {
    if (x == y)
      return true;
    else if (reachable[x].find(y) != end(reachable[x]))
      return true;
    else
      return false;
  });

  if (false) {
    cerr << "states =";
    for (auto x : states)
      cerr << " " << x;
    cerr << endl;

    cerr << "topo order =";
    for (auto x : order)
      cerr << " " << x;
    cerr << endl;
  }

  return order;
}

string HMM::get_group_consensus(size_t idx, double threshold) const {
  const string iupac = "-acmgrsvtwyhkdbn";
  string consensus = "";
  string gapped_consensus = "";

  size_t prev = 0;
  for (auto i : topological_order(transition, groups[idx].states)) {
    char present = 0;
    for (size_t j = 0; j < 4; j++)
      if (emission(i, j) >= threshold)
        present |= (1 << j);
    char cons_char = iupac[present];
    consensus += cons_char;

    if (prev != 0 and i > prev + 1)
      gapped_consensus += "[";
    if (i < prev)
      gapped_consensus += "]";
    gapped_consensus += cons_char;
    prev = i;
  }

  if (false) {
    string prev_consensus = "";
    for (auto &i : groups[idx].states) {
      char present = 0;
      for (size_t j = 0; j < 4; j++)
        if (emission(i, j) >= threshold)
          present |= (1 << j);
      prev_consensus += iupac[present];
    }

    cerr << "Consensus = " << consensus << std::endl;
    cerr << "Previous = " << prev_consensus << std::endl;
    cerr << "Gapped = " << gapped_consensus << std::endl;
  }

  return gapped_consensus;
}

string HMM::get_group_name(size_t idx) const { return groups[idx].name; }

size_t HMM::get_nstates() const { return n_states; }

size_t HMM::get_ngroups() const { return groups.size(); }

size_t HMM::get_nmotifs() const {
  size_t x = 0;
  for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if (is_motif_group(group_idx))
      x++;
  return x;
}

double HMM::get_pseudo_count() const { return pseudo_count; }

size_t HMM::get_motif_len(size_t motif_idx) const {
  size_t len = groups[motif_idx].states.size();
  return len;
}

void HMM::print_occurrence_table_header(ostream &out) const {
  out << "file\tseq\tpos\tmotifidx\tmotifname\tmotif\tstrand\tforwardpos\tcente"
         "rdist" << endl;
}

void HMM::print_occurrence_table(const string &file_path, const Data::Seq &seq,
                                 const StatePath &path, ostream &out,
                                 bool bed) const {
  size_t seqlen = seq.sequence.size();
  size_t midpoint = seqlen / 2;
  bool revcomp = false;
  // reverse-complementary sequences look like this: xxx$xxx
  // so their length is 2n + 1, and the middle nucleotide is $
  if (seqlen % 2 == 1 and seq.sequence[midpoint] == '$')
    revcomp = true;

  double center;
  if (revcomp)
    center = (midpoint - 1) / 2.0;
  else
    center = (seqlen - 1) / 2.0;

  for (size_t pos = 0; pos < path.size(); pos++)
    for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if (is_motif_group(group_idx)
          and path[pos] == groups[group_idx].states[0]) {
        size_t end = pos + 1;
        while (end != path.size() and path[end] > path[pos])
          end++;
        string motif = seq.sequence.substr(pos, end - pos);

        bool strand = (not revcomp) or (pos < midpoint);
        // forward_pos is the position relative to the forward strand
        long forward_pos = pos;
        if (strand == false)
          forward_pos = seqlen - end;
        double rel_pos = forward_pos - center;
        double motif_center_pos = rel_pos + (end - pos - 1) / 2.0;
        if (bed)
          out << seq.definition << "\t" << pos << "\t" << end << "\t"
              << groups[group_idx].name << "\t" << 0 << "\t"
              << (strand ? "+" : "-") << endl;
        else
          out << file_path << "\t" << seq.definition << "\t" << pos << "\t"
              << group_idx << "\t" << groups[group_idx].name << "\t" << motif
              << "\t" << (strand ? "+" : "-") << "\t" << forward_pos << "\t"
              << motif_center_pos << endl;
      }
}

pair<HMM, map<size_t, size_t>> HMM::add_revcomp_motifs() const {
  HMM rc = *this;
  // reverse complement motifs
  for (size_t i = 0; i < groups.size(); i++)
    if (groups[i].kind == Group::Kind::Motif)
      for (size_t j = 0; j < groups[i].states.size(); j++)
        for (size_t k = 0; k < n_emissions; k++)
          rc.emission(groups[i].states[j], k)
              = emission(groups[i].states[groups[i].states.size() - j - 1],
                         n_emissions - k - 1);

  for (size_t i = 0; i < first_state; i++) {
    rc.transition(start_state, start_state) += rc.transition(start_state, i)
        /= 2;
    rc.transition(bg_state, bg_state) += rc.transition(bg_state, i) /= 2;
  }

  HMM hmm = *this;
  hmm.add_motifs(rc);

  map<size_t, size_t> assoc;
  size_t idx = 0;
  for (size_t i = 0; i < groups.size(); i++)
    if (groups[i].kind == Group::Kind::Motif)
      assoc[i] = groups.size() + idx++;

  return make_pair(hmm, assoc);
}

namespace Exception {
namespace HMM {
namespace ParameterFile {
SyntaxError::SyntaxError(const string &token)
    : runtime_error("Syntax error: " + token) {}
UnsupportedVersion::UnsupportedVersion(size_t version)
    : runtime_error("Error: parameter file format version " + to_string(version)
                    + " not supported!") {}
}
namespace Learning {
MultipleTasks::MultipleTasks(const string &which)
    : runtime_error("Error: some " + which + " parameters are simultaneously "
                    + "assigned to multiple learning tasks.") {}
}
}
}
