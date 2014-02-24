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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <iomanip>
#include <set>
#include "../aux.hpp"
#include "hmm.hpp"

using namespace std;

double HMM::RegisteredDataSet::get_motif_prior(size_t group_idx) const {
  auto x = motif_prior.find(group_idx);
  if(x == end(motif_prior))
    throw("Error: could not find motif prior for motif group.");
  else
    return(x->second);
}

double HMM::compute_marginal_motif_prior(size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "compute_marginal_motif_prior(group_idx=" << group_idx << ")" << std::endl;
  double marginal_motif_prior = 0;
  for(auto &x: registered_datasets)
    marginal_motif_prior += x.second.class_prior * x.second.get_motif_prior(group_idx);
  return(marginal_motif_prior);
}

double HMM::get_class_motif_prior(const std::string &sha1, size_t group_idx) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "get_class_motif_prior(sha1=" << sha1 << ", group_idx=" << group_idx << ")" << std::endl;

  auto cparms = registered_datasets.find(sha1);
  if(cparms == end(registered_datasets)) {
    std::cout << "Error: failed trying to access class parameters of path with sha1: " << sha1 << "." << std::endl;
    throw("Could not find registered data set.");
  }
  auto x = cparms->second.motif_prior.find(group_idx);
  if(x == end(cparms->second.motif_prior)) {
    std::cout << "Error: failed trying to access class conditional motif prior parameters of path with sha1: " << sha1 << " for the motif with index " << group_idx << "." << std::endl;
    throw("Could not find conditional motif prior for registered data set.");
  }
  if(verbosity >= Verbosity::verbose)
    std::cout << "get_class_motif_prior(sha1=" << sha1 << ", group_idx=" << group_idx << ") = " << x->second << std::endl;
  return(x->second);
}

double HMM::get_class_prior(const std::string &sha1) const
{
  if(verbosity >= Verbosity::debug)
    std::cout << "get_class_prior(sha1=" << sha1 << ")" << std::endl;
  auto cparms = registered_datasets.find(sha1);
  if(cparms == end(registered_datasets)) {
    std::cout << "Error: failed trying to access class parameters of path with sha1: " << sha1 << "." << std::endl;
    throw("Could not find registered data set.");
  }
  if(verbosity >= Verbosity::verbose)
    std::cout << "get_class_prior(sha1=" << sha1 << ") = " << cparms->second.class_prior << std::endl;
  return(cparms->second.class_prior);
}

void HMM::shift_forward(size_t group_idx, size_t n) {
  for(size_t i = 0; i < groups[group_idx].states.size() - 1; i++)
    for(size_t j = 0; j < n_emissions; j++)
      emission(groups[group_idx].states[i],j) = emission(groups[group_idx].states[i+1],j);
  for(size_t j = 0; j < n_emissions; j++)
    emission(*groups[group_idx].states.rbegin(),j) = 1.0 / n_emissions;
  normalize_emission(emission);
}

void HMM::shift_backward(size_t group_idx, size_t n) {
  for(size_t i = groups[group_idx].states.size() - 1; i > 0; i--)
    for(size_t j = 0; j < n_emissions; j++)
      emission(groups[group_idx].states[i-1],j) = emission(groups[group_idx].states[i],j);
  for(size_t j = 0; j < n_emissions; j++)
    emission(*groups[group_idx].states.begin(),j) = 1.0 / n_emissions;
  normalize_emission(emission);
}


void HMM::serialize(ostream &os, const ExecutionInformation &exec_info, size_t format_version) const
{
  const string param_format_string = "# HMM parameter format version ";
  const string transition_matrix_string = "Transition matrix";
  const string emission_matrix_string = "Emission matrix";
  const size_t output_precision = 15;
  switch(format_version) {
    case 6:
      os << param_format_string << format_version << endl;
      os << "# " << exec_info.program_name << " " << exec_info.hmm_version << endl;
      os << "# Run on " << exec_info.datetime << endl;
      os << "# Run in " << exec_info.directory << endl;
      os << "# Command = " << exec_info.cmdline << endl;
      for(auto &x: registered_datasets) {
        os << "Dataset " << x.second.spec.path << " " << x.first << " class = " << x.second.class_prior << " motif = ";
        for(auto &y: x.second.motif_prior)
          os << " " << y.first << "/" << y.second;
        os << std::endl;
      }
      os << n_states << " states" << endl;
      os << n_emissions << " emissions" << endl;
      for(auto &motif: groups) {
        os << "Motif \"" << motif.name << "\"";
        for(auto &state: motif.states)
          os << " " << state;
        os << endl;
      }
      os << "State class";
      for(auto &motif: group_ids)
        os << " " << motif;
      os << endl;
      os << "State order";
      for(auto &o: order)
        os << " " << o;
      os << endl;
      os << transition_matrix_string << endl;
      os.precision(output_precision);
      for(size_t i = 0; i < n_states; i++) {
        os << setw(3) << fixed << i;
        for(size_t j = 0; j < n_states; j++)
          os << " " << setw(2+output_precision) << fixed << transition(i,j);
        os << endl;
      }
      os << emission_matrix_string << endl;
      for(size_t i = 0; i < n_states; i++) {
        os << setw(3) << fixed << i;
        for(size_t j = 0; j < n_emissions; j++)
          os << " " << setw(2+output_precision) << fixed << emission(i,j);
        os << endl;
      }
      break;
    default:
      throw("Format not supported!");
      break;
  }
};


void HMM::deserialize(istream &is)
{
  groups = vector<Group>();
  group_ids = vector<size_t>();
  order = vector<size_t>();

  const string param_format_string = "# HMM parameter format version ";
  string line;
  safeGetline(is, line);
  if(line.length() > param_format_string.length() and line.substr(0, param_format_string.length()) == param_format_string) {
    size_t format_version = atoi(line.substr(param_format_string.length()).c_str());
    if(format_version < 4) {
      cout << "Error loading HMM parameters: format version " << format_version << " is not supported." << endl;
    }
    while(is.good() and line[0] == '#')
      safeGetline(is,line);

    if(format_version >= 5) {
      const string key_a = "Class_prior ";
      const string key_b = "Motif_prior[0] ";
      const string key_c = "Motif_prior[1] ";
      // TODO parse information regarding registered data sets
      // if(line.substr(0, key_a.size()) == key_a)
      //   class_prior = atof(line.substr(key_a.size()).c_str());
      // else
      //   cout << "Error trying to parse HMM parameters: class prior could not be parsed: '" << line << "'." << endl;
      safeGetline(is, line);
      // TODO parse information regarding registered data sets
      // if(line.substr(0, key_b.size()) == key_b)
      //   motif_prior[0] = atof(line.substr(key_b.size()).c_str());
      // else
      //   cout << "Error trying to parse HMM parameters: first conditional motif prior could not be parsed: '" << line << "'." << endl;
      safeGetline(is, line);
      // TODO parse information regarding registered data sets
      // if(line.substr(0, key_c.size()) == key_c)
      //   motif_prior[1] = atof(line.substr(key_c.size()).c_str());
      // else
      //   cout << "Error trying to parse HMM parameters: second conditional motif prior could not be parsed: '" << line << "'." << endl;
      safeGetline(is, line);
    }

    n_states = atoi(line.c_str());
    safeGetline(is, line);
    n_emissions = atoi(line.c_str());

    if(verbosity >= Verbosity::debug)
      cout << "n_states = " << n_states << " n_emissions " << n_emissions << endl;

    safeGetline(is, line);
    while(line.substr(0,5) == "Motif") {
      size_t start = line.find('"');
      size_t stop = line.find('"', start+1);
      if(verbosity >= Verbosity::debug)
        cout << "start = " << start << " stop " << stop << endl;
      string name = line.substr(start+1, stop - start - 1);
      if(verbosity >= Verbosity::debug)
        cout << "name = '" << name << "' line " << line << endl;
      line = line.substr(stop+2);
      if(verbosity >= Verbosity::debug)
        cout << "name = " << name << " line " << line << endl;
      vector<size_t> states;
      while(line != "") {
        start = line.find(' ');
        if(verbosity >= Verbosity::debug)
          cout << "start = " << start << " line " << line << endl;
        if(start > 0) {
          size_t state = atoi(line.substr(0,start).c_str());
          states.push_back(state);
        }
        if(start == string::npos)
          break;
        else
          line = line.substr(start+1);
      }
      Group::Kind kind = Group::Kind::Motif;
      if(name == "Special")
        kind = Group::Kind::Special;
      else if(name == "Background")
        kind = Group::Kind::Background;
      Group m = {kind, name, states};
      groups.push_back(m);
      safeGetline(is, line);
    }

    // safeGetline(is, line);
    if(line.substr(0, 11) != "State class") {
      cout << "Error loading HMM parameters: expected \"State class\" but found instead: \"" << line << "\"." << endl;
      exit(-1);
    }

    line = line.substr(12);
    while(line != "") {
      size_t start = line.find(" ");
      if(start > 0) {
        size_t motif = atoi(line.substr(0,start).c_str());
        group_ids.push_back(motif);
      }
      if(start == string::npos)
        break;
      else
        line = line.substr(start+1);
    }

    safeGetline(is, line);
    if(line.substr(0, 11) != "State order") {
      cout << "Error loading HMM parameters: expected \"State order\" but found instead: \"" << line << "\"." << endl;
      exit(-1);
    }

    line = line.substr(12);
    while(line != "") {
      size_t start = line.find(" ");
      if(start > 0) {
        size_t o = atoi(line.substr(0,start).c_str());
        order.push_back(o);
      }
      if(start == string::npos)
        break;
      else
        line = line.substr(start+1);
    }


    last_state = n_states - 1;

    transition.resize(n_states, n_states);
    emission.resize(n_states, n_emissions);

    safeGetline(is,line);
    for(size_t i = 0; i < n_states; i++) {
      size_t tmp;
      is >> tmp;
      for(size_t j = 0; j < n_states; j++)
        is >> transition(i,j);
    }
    safeGetline(is,line);
    safeGetline(is,line);
    for(size_t i = 0; i < n_states; i++) {
      size_t tmp;
      is >> tmp;
      for(size_t j = 0; j < n_emissions; j++)
        is >> emission(i,j);
    }
  } else {
//    first_state = atoi(line.c_str());
    is >> transition;
    is >> emission;
    last_state = transition.size1() - 1;
  }
  if(n_states != transition.size1()) {
    cout << "Inconsistent number of states." << endl;
    exit(-1);
  }

  initialize_order_offsets();
  finalize_initialization();
};


string HMM::path2string(const HMM::StatePath &path) const
{
  string p;
  for(size_t j = 0; j < path.size(); j++)
    switch(path[j]) {
      case start_state:
        p += "^";
        break;
      default:
        if(path[j] < first_state)
          p += "0";
        else {
          char c = path[j] - first_state + 'A';
          if(c > 'Z')
            c = c - 'A' + 'a';
          if(c > 'z') {
            cout << "Error: state index is too large" << endl;
            c = '*';
          }
          p += c;
        }
        break;
    };
  return(p);
};


void HMM::to_dot(ostream &os, double minimum_transition) const
{
  static const string digits = "0123456789";
  os << "digraph d {" << endl;
  os << "  rankdir=\"LR\";" << endl;
  for(size_t i = 0; i < n_states; i++)
    for(size_t j = 0; j < n_states; j++) {
      if(transition(i,j) >= minimum_transition) {
        string label1;
        switch(i) {
          case start_state:
            label1 = "S";
            break;
          default:
            if(i < first_state)
              label1 = string("B") + digits[i-bg_state];
            else
              label1 = string() + digits[i-first_state];
            break;
        }
        string label2;
        switch(j) {
          case start_state:
            label2 = "S";
            break;
          default:
            if(j < first_state)
              label2 = string("B") + digits[j-bg_state];
            else
              label2 = string() + digits[j-first_state];
            break;
        }
        os << label1 << " -> " << label2 << ";" << endl;
      }
    }
  os << "}" << endl;
};


bool HMM::check_consistency_transitions() const
{
  for(size_t k = 0; k < n_states; k++) {
    double p = 0;
    for(size_t l = 0; l < n_states; l++)
      p += transition(k,l);
    if(fabs(1-p)>1e-6)
      return(false);
  }
  return(true);
}


bool HMM::check_consistency_emissions() const
{
  // TODO: implement this correctly w.r.t. higher order emissions
  /*
  for(size_t k = 0; k < n_states; k++) {
    double p = 0;
    for(size_t l = 0; l < n_emissions; l++)
      p += emission(k,l);
    if(fabs(1-p)>1e-6)
      return(false);
  }  */
  return(true);
}

bool HMM::check_consistency() const
{
  bool ok = true;
  ok = ok and check_consistency_transitions();
  ok = ok and check_consistency_emissions();
  ok = ok and (group_ids.size() == order.size() and group_ids.size() == n_states);
  for(auto &m: groups)
    for(auto &state: m.states)
      if(state > n_states)
        ok = false;
  for(auto &o: order)
    if(o > max_order)
      ok = false;
  for(auto &d: group_ids)
    if(d > groups.size())
      ok = false;
  return(ok);
}

ostream &operator<<(ostream& os, const HMM &hmm)
{
  os << "HMM: " << hmm.n_states << " states (including start state), " << hmm.n_emissions << " emissions." << endl
    << "Transition probability matrix: " << hmm.transition << endl
    << "Emission probability matrix: " << hmm.emission << endl;
  return(os);
}


/** Set the emissions according to the desired matrix.
 * Also embed the matrix with pad_left uniform emission positions on the left,
 * and pad_left uniform emission positions on the right. Finally, a number of
 * n_insertions positions with uniform emissions are added for insertions states.
 */
void HMM::set_motif_emissions(const matrix_t &e, size_t first_padded, size_t n_insertions, size_t pad_left, size_t pad_right)
{
  if(verbosity >= Verbosity::debug)
    cout << "set motif_emissions; n_insertions = " << n_insertions << endl;

  const size_t n = e.size1();

  // first set the emissions of the padding states on the left to uniform
  for(size_t i = 0; i < pad_left; i++) {
    for(size_t j = 0; j < e.size2(); j++)
      emission(first_padded + i , j) = 1.0 / e.size2();
    for(size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + i , j) = 0;
  }

  // then set the emissions of the states of the actual motif to the given matrix
  for(size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + i , j) = e(i, j);
    for(size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + i , j) = 0;
  }

  // finally set the emissions of the padding states on the right to uniform
  for(size_t i = 0; i < pad_right; i++) {
    for(size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + n + i , j) = 1.0 / e.size2();
    for(size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + n + i , j) = 0;
  }

  // set the emissions of the insert states to uniform
  for(size_t i = 0; i < n_insertions; i++) {
    for(size_t j = 0; j < e.size2(); j++)
      emission(first_padded + pad_left + n + pad_right + i, j) = 1.0 / e.size2();
    for(size_t j = e.size2(); j < n_emissions; j++)
      emission(first_padded + pad_left + n + pad_right + i, j) = 0;
  }
}

// NOTE: this assumes that the alphabet has size 4, which is okay for nucleic acids, but not for amino acids
void HMM::update_history(History &history, size_t obs) const
{
  if(verbosity >= Verbosity::everything)
    cout << "obs = " << obs << endl
      << "history = " << history.observation << endl
      << "history length = " << history.length << endl
      << "alphabet_size = " << alphabet_size << endl
      << "max_order  = " << max_order  << endl
      << "history * alphabet_size = " << history.observation * alphabet_size << endl
      << "history << 2 = " << (history.observation << 2) << endl
      << "ipow(alphabet_size, max_order + 1) = " << ipow(alphabet_size, max_order + 1) << endl
      << "(1 << (2 * (max_order + 1))) = " << (1 << (2 * (max_order + 1))) << endl;
  if(obs == empty_symbol) {
    history.observation = 0;
    history.length = 0;
  } else {
    history.observation = ((history.observation << 2)         // left shift one position
        + obs)                                                // add the new observation
      % (1 << (2 * (max_order + 1)));                         // truncate to max_order position
    history.length++;
    history.length = min<size_t>(history.length, max_order + 1);
  }
  if(verbosity >= Verbosity::everything)
    cout << "new history = " << history.observation << " length = " << history.length << endl;
}

// NOTE: this assumes that the alphabet has size 4, which is okay for nucleic acids, but not for amino acids
void HMM::update_history_front(History &history, size_t obs) const
{
  if(verbosity >= Verbosity::everything)
    cout << "adding obs = " << obs << endl
      << "prev history = " << history.observation << endl
      << "prev history length = " << history.length<< endl
      << "alphabet_size = " << alphabet_size << endl
      << "max_order  = " << max_order  << endl
      << "history / alphabet_size = " << history.observation / alphabet_size << endl
      << "history >> 2 = " << (history.observation >> 2) << endl;
  if(obs == empty_symbol) {
    if(history.length > 0) {
      history.observation = history.observation >> 2;
      history.length--;
    }
  } else {
    history.observation = (obs << (2 * max_order))    // left shift obs to the most significant position
      + (history.observation >> 2);                   // add the previous history, right-shifted one position
    history.length++;
    history.length = min<size_t>(history.length, max_order + 1);
  }
  if(verbosity >= Verbosity::everything)
    cout << "new history = " << history.observation << endl;
}

void HMM::initialize_order_offsets()
{
  max_order = 0;
  for(auto &o: order)
    if(o > max_order)
      max_order = o;

  order_offset.resize(max_order);
  double n = 0;
  for(size_t i = 0; i < max_order; i++)
    order_offset[i] = n += ipow(alphabet_size, max_order - i + 1);
  if(verbosity >= Verbosity::debug) {
    cout << "order =";
    for(auto o: order)
      cout << " " << o;
    cout << endl;
    cout << "order_offset =";
    for(auto o: order_offset)
      cout << " " << o;
    cout << endl;
  }
}

void HMM::normalize_transition(matrix_t &m) const
{
  for(size_t k = 0; k < m.size1(); k++) {
    double z = 0;
    for(size_t l = 0; l < m.size2(); l++)
      z += m(k,l);
    if(z > 0)
      for(size_t l = 0; l < m.size2(); l++)
        m(k,l) = m(k,l) / z;
    else
      for(size_t l = 0; l < m.size2(); l++)
        m(k,l) = 0;
  }
}

void HMM::normalize_emission(matrix_t &m) const
{
  for(size_t k = bg_state; k < m.size1(); k++) {
    size_t offset = 0;
    if(verbosity >= Verbosity::debug)
      cout << "k = " << k << endl
        << "order[" << k << "] = " << order[k] << endl;
    for(int i = order[k]; i >= 0; i--) {
      size_t begin = offset;
      offset += ipow(alphabet_size, i+1);
      size_t end = offset;

      if(verbosity >= Verbosity::debug)
        cout << "Doing all emissions in " << begin << " to " << end << endl;
      size_t x = begin;
      while(x < end) {
        size_t y = x + alphabet_size;
        if(verbosity >= Verbosity::debug)
          cout << "Doing emission " << x << " to " << y << endl;

        double z = 0;
        for(size_t b = x; b < y; b++)
          z += m(k,b);
        if(z > 0)
          for(size_t b = x; b < y; b++)
            m(k,b) = m(k,b) / z;
        else
          for(size_t b = x; b < y; b++)
            m(k,b) = 0;

        x = y;
      }
    }
  }

  // enforce emissions of start state
  for(size_t b = 0; b < n_emissions; b++)
    m(start_state,b) = 0;
}

Training::Tasks HMM::define_training_tasks(const hmm_options &options) const
{
  Training::Tasks tasks;

  // First the (mostly) discriminative tasks
  bool atleast_one_discriminative_task = false;
  bool do_transition = options.bg_learning == Training::Method::Gradient and options.objectives.size() == 1;
  for(auto &objective: options.objectives) {
    if(not Measures::is_discriminative(objective.measure))
      continue;
    atleast_one_discriminative_task = true;
    Training::Task task;
    task.motif_name = objective.motif_name;
    task.measure = objective.measure;
    task.series_expression = objective.series_expression;
    for(auto &group: groups) {
      if(group.name == objective.motif_name) {
        copy(begin(group.states), end(group.states), back_inserter(task.targets.emission));
        if(do_transition)
          copy(begin(group.states), end(group.states), back_inserter(task.targets.transition));
      }
    }
    tasks.push_back(task);

    if(options.verbosity >= Verbosity::verbose) {
      cout << "Generated primary training targets." << endl << "Emissions:";
      for(auto e: task.targets.emission)
        cout << " " << e;
      cout << endl;
      cout << "Transitions:";
      for(auto t: task.targets.transition)
        cout << " " << t;
      cout << endl;
    }
  }

  // get all assigned emission and transition states, and check that none are doubly assigned
  set<size_t> e_states, t_states;
  for(auto &task: tasks) {
    for(auto &e_state: task.targets.emission) {
      auto pair = e_states.insert(e_state);
      if(not pair.second) {
        cout << "Error: some emission parameters are supposed to be simultaneously subjected to multiple learning tasks." << endl;
        exit(-1);
      }
    }
    for(auto &t_state: task.targets.transition) {
      auto pair = t_states.insert(t_state);
      if(not pair.second) {
        cout << "Error: some transition parameters are supposed to be simultaneously subjected to multiple learning tasks." << endl;
        exit(-1);
      }
    }
  }

  // get all series names
  set<string> series_names;
  for(auto &task: tasks)
    for(auto &expr: task)
      series_names.insert(expr.series);

  // And then the generative part
  if(options.bg_learning != Training::Method::None) {
    Training::Task task;
    if(atleast_one_discriminative_task) {
      task.motif_name = "Generative parameters";
      task.measure = Measure::Likelihood;
    } else {
      task.motif_name = options.objectives[0].motif_name;
      task.measure = options.objectives[0].measure;
    }

    for(auto &series_name: series_names)
      task.series_expression.push_back({+1, series_name});
    for(size_t i = 0; i < n_states; i++) {
      if(e_states.find(i) == end(e_states))
        if(i > start_state)
          task.targets.emission.push_back(i);
      if(t_states.find(i) == end(t_states))
        task.targets.transition.push_back(i);
    }

    tasks.push_back(task);

    if(options.verbosity >= Verbosity::verbose) {
      cout << "Generated secondary training targets." << endl << "Emissions:";
      for(auto e: task.targets.emission)
        cout << " " << e;
      cout << endl;
      cout << "Transitions:";
      for(auto t: task.targets.transition)
        cout << " " << t;
      cout << endl;
    }
  }
  return(tasks);
}

size_t HMM::non_zero_parameters(const Training::Targets &training_targets) const
{
  size_t z = 0;
  for(auto i: training_targets.transition)
    for(size_t j = 0; j < n_states; j++)
      if(transition(i,j) > 0)
        z++;
  for(auto i: training_targets.emission)
    for(size_t j = 0; j < n_emissions; j++)
      if(emission(i,j) > 0)
        z++;
  if(verbosity >= Verbosity::verbose)
    cout << "non_zero_parameters = " << z << endl;
  return(z);
}

size_t HMM::n_parameters() const
{
  size_t n = (n_states * (n_states - 1)) + (n_states * (3 * ipow(alphabet_size, max_order)));
  if(verbosity >= Verbosity::verbose)
    cout << "n_parameters = " << n << endl;
  return(n);
}

size_t HMM::count_motif(const HMM::StatePath &path, size_t motif) const
{
  if(motif >= groups.size())
    return(0);
  size_t first_state = *groups[motif].states.begin(); // TODO make this more flexible
  size_t n = 0;
  for(auto x: path)
    if(x == first_state)
      n++;
  return(n);
}

bool HMM::is_motif_group(size_t group_idx) const
{
  return(groups[group_idx].kind == Group::Kind::Motif);
}

bool HMM::is_motif_state(size_t state_idx) const
{
  return(is_motif_group(group_ids[state_idx]));
}

HMM::mask_t HMM::compute_mask(const Data::Collection &data) const
{
  if(verbosity >= Verbosity::debug)
    cout << "HMM::compute-mask(Data::Series)" << endl;
  HMM::mask_t mask;
  for(auto &series: data)
    for(auto &data_set: series) {
      mask_sub_t m;
      for(auto &seq: data_set) {
        vector<size_t> v;
        HMM::StatePath path;
        viterbi(seq, path);
        size_t idx = 0;
        for(auto state: path) {
          if(state >= first_state)
            v.push_back(idx); // TODO: check: should I insert idx - 1 instead?
                              // -> No; but there was an issue with reverse complements; fixed in src/plasma/data.cpp DataSet::mask by pos = 2 * n - pos where n is the sequence length

          idx++;
        }
        if(not v.empty())
          m[seq.definition] = v;
      }
      if(not m.empty())
        mask[data_set.path] = m;
    }
  return(mask);
}

Training::Range HMM::complementary_states(size_t group_idx) const
{
  Training::Range range;
  for(size_t i = start_state; i < n_states; i++)
    if(not is_motif_state(i) or group_ids[i] != group_idx)
      range.push_back(i);
  if(verbosity >= Verbosity::debug) {
    cout << "Complementary states =";
    for(auto &x: range)
      cout << " " << x;
    cout << endl;
  }
  return(range);
}

string HMM::get_group_consensus(size_t idx, double threshold) const
{
  const string iupac = "-acmgrsvtwyhkdbn";
  string consensus = "";
  for(auto &i: groups[idx].states) {
    char present = 0;
    for(size_t j = 0; j < 4; j++)
      if(emission(i,j) >= threshold)
        present |= (1 << j);
    consensus += iupac[present];
  }
  return(consensus);
}

string HMM::get_group_name(size_t idx) const
{
  return(groups[idx].name);
}

size_t HMM::get_ngroups() const
{
  return(groups.size());
}

size_t HMM::get_nmotifs() const
{
  size_t x = 0;
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if(is_motif_group(group_idx))
      x++;
  return(x);
}

double HMM::get_pseudo_count() const
{
  return(pseudo_count);
}

size_t HMM::get_motif_len(size_t motif_idx) const
{
  size_t len = groups[motif_idx].states.size();
  return(len);
}

