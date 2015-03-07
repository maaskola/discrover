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
 *       Filename:  hmm_learn.cpp
 *
 *    Description:  Learning routines for the HMM class
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <fstream>
#include <iomanip>
#include "../timer.hpp"
#include "../aux.hpp"
#include "hmm.hpp"
#include "../format_constants.hpp"

using namespace std;

#define DO_PARALLEL 1

string line_search_status(int status) {
  string msg;
  switch (status) {
    case 0:
      msg = "IMPROPER INPUT PARAMETERS.";
      break;
    case 1:
      msg = "THE SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION HOLD.";
      break;
    case -1:
      msg = "A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.";
      break;
    case 2:
      msg = "RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.";
      break;
    case 3:
      msg = "NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.";
      break;
    case 4:
      msg = "THE STEP IS AT THE LOWER BOUND STPMIN.";
      break;
    case 5:
      msg = "THE STEP IS AT THE UPPER BOUND STPMAX.";
      break;
    case 6:
      msg = "ROUNDING ERRORS PREVENT FURTHER PROGRESS. THERE MAY NOT BE A STEP WHICH SATISFIES THE SUFFICIENT DECREASE AND CURVATURE CONDITIONS. TOLERANCES MAY BE TOO SMALL.";
      break;
  }
  return msg;
}

void HMM::train_background(const Data::Collection &collection,
                           const Options::HMM &options) {
  HMM bg_hmm(options.verbosity);

  Training::Task task;
  task.measure = Measure::Likelihood;

  for (size_t i = bg_hmm.start_state; i < bg_hmm.n_states; i++) {
    task.targets.transition.push_back(i);
    if (i != start_state)
      task.targets.emission.push_back(i);
  }

  bg_hmm.reestimation(collection, task, options);

  if (options.verbosity >= Verbosity::debug)
    cout << "Resulting model = " << bg_hmm << endl << endl;

  for (size_t i = 0; i < bg_hmm.transition.size1(); i++)
    for (size_t j = 0; j < bg_hmm.transition.size2(); j++)
      transition(i, j) = bg_hmm.transition(i, j);
  for (size_t i = 0; i < bg_hmm.emission.size1(); i++)
    for (size_t j = 0; j < bg_hmm.emission.size2(); j++)
      emission(i, j) = bg_hmm.emission(i, j);

  if (transition.size1() > bg_hmm.transition.size1())
    throw Exception::HMM::Learning::TrainBgTooLate();

  if (options.verbosity >= Verbosity::debug)
    cout << *this << endl;
}

void HMM::initialize_bg_with_bw(const Data::Collection &collection,
                                const Options::HMM &options) {
  if (options.verbosity >= Verbosity::info)
    cout << "Initializing background with Baum-Welch algorithm." << endl;

  Options::HMM bg_options = options;
  if (options.verbosity == Verbosity::info)
    bg_options.verbosity = Verbosity::error;

  Timer timer;
  train_background(collection, bg_options);
  double time = timer.tock();

  if (options.timing_information)
    cerr << "Background learning: " + time_to_pretty_string(time) << endl;
}

Training::Result HMM::train(const Data::Collection &collection,
                            const Training::Tasks &tasks,
                            const Options::HMM &options) {
  Training::Result result;
  if (options.verbosity >= Verbosity::verbose)
    cout << "Model to be evaluated = " << *this << endl;
  if (tasks.empty()) {
    if (options.verbosity >= Verbosity::info)
      cout
          << "Not performing training because no training tasks were specified."
          << endl;
  } else {
    bool any_found = false;
    for (auto &group : groups)
      for (auto &task : tasks)
        if (task.motif_name == group.name) {
          any_found = true;
        }
    if (not any_found) {
      if (options.verbosity >= Verbosity::info)
        cout << "Skipping training because no motifs specified in the tasks "
                "have corresponding states in the HMM." << endl;
    } else {
      if (options.verbosity >= Verbosity::info)
        cout << "Performing training." << endl;
      Timer learning_timer;

      if (options.verbosity >= Verbosity::verbose)
        cout << "Registering data sets for class based HMMs." << endl;

      for (auto &contrast : collection)
        for (auto &dataset : contrast)
          registration.add_dataset(
              dataset, (1.0 * dataset.set_size) / collection.set_size);
      for (auto &task : tasks)
        registration.add_bitmask(task.motif_name, compute_bitmask(task),
                                 options.conditional_motif_prior1,
                                 options.conditional_motif_prior2);

      result = train_inner(collection, tasks, options);
      if (options.verbosity >= Verbosity::verbose)
        cout << endl << "The parameters changed by an L1-norm of "
             << result.delta << endl;

      double time = learning_timer.tock();
      if (options.timing_information)
        cerr << "Learning: " + time_to_pretty_string(time) << endl;

      if (options.verbosity >= Verbosity::debug)
        cout << "HMM after training:" << endl << *this << endl;

      result.parameter_file = options.label + ".hmm";

      if (options.verbosity >= Verbosity::info)
        cout << endl << left << setw(report_col_width) << "HMM parameters"
             << right << result.parameter_file << endl;

      ofstream os(result.parameter_file.c_str());
      serialize(os, options.exec_info);
    }
  }
  return result;
}

Training::Result HMM::train_inner(const Data::Collection &collection,
                                  const Training::Tasks &tasks,
                                  const Options::HMM &options) {
  Training::Result result;
  if (tasks.empty())
    return result;
  else {
    HMM previous(*this);

    if (options.sampling.do_sampling) {
      Training::Task task
          = tasks[0];  // TODO make HMM sampling work with multiple tasks
      auto sampling_results = mcmc(collection, task, options);
      if (options.verbosity >= Verbosity::info)
        cout << "Sampling done. Got " << sampling_results.size() << " results."
             << endl;
      double max = -numeric_limits<double>::infinity();
      for (auto y : sampling_results)
        for (auto x : y) {
          if (options.verbosity >= Verbosity::verbose)
            cout << "Checking if we have a new maximum: " << x.second << endl;
          if (x.second > max) {
            *this = x.first;
            max = x.second;
            if (options.verbosity >= Verbosity::info)
              cout << "Found new maximum: " << max << endl;
          }
        }
    } else
      result.state = iterative_training(collection, tasks, options);
    result.delta = norml1(previous.emission - emission)
                   + norml1(previous.transition - transition);
    return result;
  }
}

double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E,
                               const Data::Collection &collection,
                               const Training::Targets &targets,
                               const Options::HMM &options) const {
  double log_likel = 0;
  for (auto &contrast : collection)
    for (auto &dataset : contrast)
      if (not dataset.is_control)
        log_likel += BaumWelchIteration(T, E, dataset, targets, options);
  if (verbosity >= Verbosity::debug)
    cout << "Done BaumWelchIteration(Collection) log_likel = " << log_likel
         << endl;
  return log_likel;
}

double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E,
                               const Data::Set &dataset,
                               const Training::Targets &targets,
                               const Options::HMM &options) const {
  double log_likel = 0;
  // TODO find out if there's a performance benefit to the split-up critical
  // regions.
#pragma omp parallel for reduction(+ : log_likel) shared(E, T) if (DO_PARALLEL)
  for (size_t j = 0; j < dataset.sequences.size(); j++)
    log_likel += BaumWelchIteration(T, E, dataset.sequences[j], targets);
  if (verbosity >= Verbosity::debug)
    cout << "Done BaumWelchIteration(Seqs) log_likel = " << log_likel << endl;
  return log_likel;
}

double HMM::BaumWelchIteration_single(matrix_t &T, matrix_t &E,
                                      const Data::Seq &s,
                                      const Training::Targets &targets) const {
  size_t L = s.isequence.size();

  vector_t scale;
  matrix_t f = compute_forward_scaled(s, scale);
  matrix_t b = compute_backward_prescaled(s, scale);

  double log_likel = log_likelihood_from_scale(scale);

  if (not targets.transition.empty()) {
    if (not(T.size1() == n_states and T.size2() == n_states))
      T = zero_matrix(n_states, n_states);

    // for all transitions except the one to the start state
    for (size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if (symbol == empty_symbol)
        for (auto k : targets.transition)
          T(k, start_state) += f(i, k) * transition(k, start_state)
                               * b(i + 1, start_state);
      else
        // TODO change the order of the loops
        for (auto k : targets.transition) {
          double f_i_k = f(i, k);
          for (auto suc : succ[k])
            T(k, suc) += f_i_k * transition(k, suc) * emission(suc, symbol)
                         * b(i + 1, suc);
        }
    }

    // for the transition to the start_state
    for (auto pre : pred[start_state])
      T(pre, start_state) += f(L, pre) * transition(pre, start_state)
                             * b(L + 1, start_state);
  }

  if (not targets.emission.empty()) {
    if (not(E.size1() == n_states and E.size2() == n_emissions))
      E = zero_matrix(n_states, n_emissions);

    for (size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if (symbol != empty_symbol)
        for (auto k : targets.emission)
          E(k, symbol) += f(i + 1, k) * b(i + 1, k) * scale(i + 1);
    }
  }
  return log_likel;
}

double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Seq &s,
                               const Training::Targets &targets) const {
  size_t L = s.isequence.size();

  vector_t scale;
  matrix_t f = compute_forward_scaled(s, scale);
  matrix_t b = compute_backward_prescaled(s, scale);

  double log_likel = log_likelihood_from_scale(scale);
  if (verbosity >= Verbosity::debug)
    cerr << "log_likel = " << log_likel << endl;

  if (not targets.transition.empty()) {
    matrix_t t = zero_matrix(n_states, n_states);

    // for all transitions except the one to the start state
    for (size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if (symbol == empty_symbol)
        for (auto k : targets.transition)
          t(k, start_state) += f(i, k) * transition(k, start_state)
                               * b(i + 1, start_state);
      else
        // TODO change the order of the loops
        for (auto k : targets.transition) {
          double f_i_k = f(i, k);
          for (auto suc : succ[k])
            t(k, suc) += f_i_k * transition(k, suc) * emission(suc, symbol)
                         * b(i + 1, suc);
        }
    }

    // for the transition to the start_state
    for (auto pre : pred[start_state])
      t(pre, start_state) += f(L, pre) * transition(pre, start_state)
                             * b(L + 1, start_state);

    if (verbosity >= Verbosity::debug)
      cerr << "t = " << t << endl;

#pragma omp critical(update_T)
    T += t;
  }

  if (not targets.emission.empty()) {
    matrix_t e = zero_matrix(n_states, n_emissions);

    for (size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if (symbol != empty_symbol)
        for (auto k : targets.emission)
          e(k, symbol) += f(i + 1, k) * b(i + 1, k) * scale(i + 1);
    }

    if (verbosity >= Verbosity::debug)
      cerr << "e = " << e << endl;

#pragma omp critical(update_E)
    E += e;
  }
  if (verbosity >= Verbosity::debug)
    cout << "Done BaumWelchIteration(Seq) log_likel = " << log_likel << endl;
  return log_likel;
}

double HMM::ViterbiIteration(matrix_t &T, matrix_t &E,
                             const Data::Collection &collection,
                             const Training::Targets &targets,
                             const Options::HMM &options) {
  double log_likel = 0;
  for (auto &contrast : collection)
    for (auto &dataset : contrast)
      log_likel += ViterbiIteration(T, E, dataset, targets, options);
  return log_likel;
}

double HMM::ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Set &dataset,
                             const Training::Targets &training_targets,
                             const Options::HMM &options) {
  double log_likel = 0;
#pragma omp parallel for reduction(+ : log_likel) shared(E, T) if (DO_PARALLEL)
  for (size_t j = 0; j < dataset.sequences.size(); j++) {
    StatePath path;
    double cur_log_likel = viterbi(dataset.sequences[j], path);

    size_t L = dataset.sequences[j].isequence.size();

    matrix_t t = zero_matrix(n_states, n_states);
    if (not training_targets.transition.empty()) {
      t(start_state, path[0]) += 1;
      for (size_t i = 0; i < L - 1; i++)
        t(path[i], path[i + 1]) += 1;
      t(path[L - 1], start_state) += 1;
    }

    matrix_t e = zero_matrix(n_states, n_emissions);
    if (not training_targets.emission.empty())
      for (size_t i = 0; i < L; i++)
        e(path[i], dataset.sequences[j].isequence[i]) += 1;

#pragma omp critical
    {
      if (not training_targets.emission.empty())
        E += e;
      if (not training_targets.transition.empty())
        T += t;
    }
    log_likel += cur_log_likel;
  }
  return log_likel;
}

void HMM::reestimation(const Data::Collection &collection,
                       const Training::Task &task,
                       const Options::HMM &options) {
  Training::Tasks tasks;
  tasks.push_back(task);
  iterative_training(collection, tasks, options);
}

double HMM::reestimationIteration(const Data::Collection &collection,
                                  const Training::Task &task,
                                  const Options::HMM &options) {
  double log_likel = -numeric_limits<double>::infinity();
  if (store_intermediate)
    serialize(cerr, options.exec_info);
  matrix_t T, E;
  if (not task.targets.transition.empty())
    T = zero_matrix(n_states, n_states);
  if (not task.targets.emission.empty())
    E = zero_matrix(n_states, n_emissions);

  switch (task.measure) {
    case Measure::Likelihood:
      log_likel = BaumWelchIteration(T, E, collection, task.targets, options);
      break;
    case Measure::Viterbi:
      log_likel = ViterbiIteration(T, E, collection, task.targets, options);
      break;
    default:
      break;
  }

  M_Step(T, E, task.targets, options);

  if (verbosity >= Verbosity::verbose)
    cout << "Done HMM::reestimationIteration log_likel = " << log_likel << endl;
  return log_likel;
};

void HMM::M_Step(const matrix_t &T_, const matrix_t &E_,
                 const Training::Targets &targets,
                 const Options::HMM &options) {
  if (verbosity >= Verbosity::verbose) {
    if (not targets.transition.empty())
      cerr << "Before normalization: New transitions = " << T_ << endl;
    if (not targets.emission.empty())
      cerr << "Before normalization: New emissions = " << E_ << endl;
  }

  matrix_t T, E;
  // add pseudo count to expected transitions
  if (not targets.transition.empty()) {
    T = T_;
    for (auto t : targets.transition)
      for (size_t j = 0; j < n_states; j++)
        if (transition(t, j) > 0)
          T(t, j) += options.transition_pseudo_count;
  }

  // add pseudo count to expected emissions
  if (not targets.emission.empty()) {
    E = E_;
    for (auto t : targets.emission)
      for (size_t j = 0; j < n_emissions; j++)
        E(t, j) += options.emission_pseudo_count;
  }

  // normalize transition probabilities
  if (not targets.transition.empty())
    normalize_transition(T);

  // normalize emission probabilities of all but the start state
  if (not targets.emission.empty())
    normalize_emission(E);

  for (auto t : targets.transition)
    for (size_t j = 0; j < n_states; j++)
      transition(t, j) = T(t, j);
  for (auto t : targets.emission)
    for (size_t j = 0; j < n_emissions; j++)
      emission(t, j) = E(t, j);

  if (verbosity >= Verbosity::verbose) {
    if (not targets.transition.empty())
      cerr << "New transitions = " << transition << endl;
    if (not targets.emission.empty())
      cerr << "New emissions = " << emission << endl;
  }
}

Gradient HMM::compute_gradient(const Data::Collection &collection,
                               double &score, const Training::Task &task,
                               bool weighting) const {
  if (verbosity >= Verbosity::verbose) {
    cerr << "HMM::compute_gradient(Data::Collection)" << endl
         << "Task = " << task.motif_name << ":";
    for (auto &expr : task.contrast_expression)
      cerr << expr;
    cerr << ":" << task.measure << endl;
  }
  Gradient gradient;
  if (not task.targets.transition.empty() and gradient.transition.size1() == 0)
    gradient.transition = zero_matrix(n_states, n_states);
  if (not task.targets.emission.empty() and gradient.emission.size1() == 0)
    gradient.emission = zero_matrix(n_states, n_emissions);

  score = 0;
  double W = 0;
  for (auto &expr : task.contrast_expression) {
    auto iter = collection.find(expr.contrast);
    if (iter != end(collection)) {
      Gradient g;
      double s = compute_gradient(*iter, g, task);
      double w = iter->set_size;
      W += w;
      if (not task.targets.transition.empty())
        gradient.transition += g.transition * expr.sign * (weighting ? w : 1.0);
      if (not task.targets.emission.empty())
        gradient.emission += g.emission * expr.sign * (weighting ? w : 1.0);
      if (verbosity >= Verbosity::verbose)
        cerr << "contrast = " << expr.contrast << " score = " << s << endl;
      score += expr.sign * s * (weighting ? w : 1.0);
    }
  }
  if (weighting) {
    score /= W;
    if (not task.targets.transition.empty())
      gradient.transition /= W;
    if (not task.targets.emission.empty())
      gradient.emission /= W;
  }
  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::compute_gradient::end score = " << score << endl;
  return gradient;
}

Gradient HMM::compute_gradient(const Data::Contrast &contrast, double &score,
                               const Training::Task &task) const {
  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::compute_gradient(Data::Contrast, task)" << endl;
  // ObjectiveFunction objFunc = objFunc_for_training(task.method);
  Gradient gradient;
  score = 0;
  for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if (is_motif_group(group_idx))
      score += compute_gradient(contrast, gradient, task);
  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::compute_gradient(Data::Contrast, task)::end" << endl;
  return gradient;
}

double HMM::compute_gradient(const Data::Contrast &contrast, Gradient &gradient,
                             const Training::Task &task) const {
  bitmask_t present = compute_bitmask(task);

  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::compute_gradient(Data::Contrast, task, present)" << endl;
  double score = 0;
  switch (task.measure) {
    case Measure::MutualInformation:
      score = mutual_information_gradient(contrast, task, present, gradient);
      break;
    case Measure::RankInformation:
      score = rank_information_gradient(contrast, task, present, gradient);
      break;
    case Measure::ChiSquare:
      score = chi_square_gradient(contrast, task, present, gradient);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      score = matthews_correlation_coefficient_gradient(contrast, task, present,
                                                        gradient);
      break;
    case Measure::LogLikelihoodDifference:
      score = log_likelihood_difference_gradient(contrast, task, present,
                                                 gradient);
      break;
    case Measure::DeltaFrequency:
      score = site_frequency_difference_gradient(contrast, task, present,
                                                 gradient);
      break;
    case Measure::ClassificationPosterior:
    case Measure::ClassificationLikelihood:
      score = class_likelihood_gradient(contrast, task, present, gradient);
      break;
    default:
      throw Exception::HMM::Learning::GradientNotImplemented(task.measure);
  }

  if (verbosity >= Verbosity::verbose)
    cout << "Gradient calculation yielded a score of " << score << "." << endl;
  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::compute_gradient(Data::Contrast, task, present)::end" << endl;

  return score;
}

Training::State HMM::iterative_training(const Data::Collection &collection,
                                        const Training::Tasks &tasks,
                                        const Options::HMM &options) {
  Training::State state(tasks.size());
  size_t iteration = 0;

  if (verbosity >= Verbosity::info) {
    cout << endl << "Iteration                                      "
         << iteration << endl;
    for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if (is_motif_group(group_idx))
        cout << "Motif consensus                                "
             << groups[group_idx].name << ":" << get_group_consensus(group_idx)
             << endl;
  }

  while ((iteration++ < options.termination.max_iter
          or options.termination.max_iter == 0)
         and perform_training_iteration(collection, tasks, options, state))
    if (verbosity >= Verbosity::info) {
      cout << endl << "Iteration                                      "
           << iteration << endl;
      // cout << "Gradient learning, relative score difference   " <<
      // relative_score_difference << endl;
      for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
        if (is_motif_group(group_idx))
          cout << "Motif consensus                                "
               << groups[group_idx].name << ":"
               << get_group_consensus(group_idx) << endl;
    }

  if (verbosity >= Verbosity::info) {
    cout << endl << "Finished after " << iteration << " iterations." << endl
         << "Final score" << (tasks.size() > 1 ? "s" : "") << " = ";
    for (size_t i = 0; i < tasks.size(); i++)
      cout << " " << to_string(tasks[i]) << "=" << *state.scores[i].rbegin();
    cout << endl;
    for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if (is_motif_group(group_idx))
        cout << " " << groups[group_idx].name << ":"
             << get_group_consensus(group_idx);
    cout << endl;
  }

  return state;
}

bool HMM::perform_training_iteration(const Data::Collection &collection,
                                     const Training::Tasks &tasks,
                                     const Options::HMM &options,
                                     Training::State &state) {
  bool done = true;

  size_t task_idx = 0;
  for (auto task : tasks) {
    double score = -numeric_limits<double>::infinity();
    if (task.measure != Measure::Undefined) {
      if (store_intermediate)
        serialize(cerr, options.exec_info);

      if (verbosity >= Verbosity::verbose)
        cerr << "Current parameters: " << *this << endl;

      // collect the states whose parameters are subject to this task
      Training::Range task_targets(task.targets.transition.size()
                                   + task.targets.emission.size());
      copy(begin(task.targets.transition), end(task.targets.transition),
           begin(task_targets));
      copy(begin(task.targets.emission), end(task.targets.emission),
           begin(task_targets));
      sort(begin(task_targets), end(task_targets));
      task_targets.resize(distance(
          begin(task_targets), unique(begin(task_targets), end(task_targets))));

      // check if the task states are reachable
      if (find(begin(task_targets), end(task_targets), start_state)
          == end(task_targets)) {
        // the task states do not include the start state
        // check if the task states are reachable from the non-task states
        double z = 0;
        for (size_t i = 0; i < n_states; i++)
          if (find(begin(task_targets), end(task_targets), i)
              == end(task_targets))
            for (auto j : task_targets)
              z += transition(i, j);
        if (z == 0) {
          if (verbosity >= Verbosity::info)
            cout << "None of the task states are reachable from the non-task "
                    "states; skipping task." << endl;
          continue;
        }
      }

      if (options.termination.past <= state.scores[task_idx].size())
        score
            = *(state.scores[task_idx].rbegin() + options.termination.past - 1);

      if (Training::measure2method(task.measure) == Training::Method::Gradient)
        done = perform_training_iteration_gradient(
                   collection, task, options, state.center, score) and done;

      if (Training::measure2method(task.measure)
          == Training::Method::Reestimation)
        done = perform_training_iteration_reestimation(collection, task,
                                                       options, score) and done;

      if ((task.measure == Measure::ClassificationPosterior
           or task.measure == Measure::ClassificationLikelihood)
          and (not(options.dont_learn_class_prior
                   and options.dont_learn_conditional_motif_prior)))
        done = reestimate_class_parameters(collection, task, options, score)
               and done;
    }
    state.scores[task_idx++].push_back(score);
  }

  if (store_intermediate)
    serialize(cerr, options.exec_info);

  return not done;
}

bool HMM::reestimate_class_parameters(const Data::Collection &collection,
                                      const Training::Task &task,
                                      const Options::HMM &options,
                                      double &score) {
  if (verbosity >= Verbosity::debug)
    cout << "HMM::reestimate_class_parameters(Data::Collection, ...)" << endl;

  bitmask_t present = compute_bitmask(task);

  unordered_map<string, double> class_counts;
  unordered_map<string, double> motif_counts;
  for (auto &contrast : collection)
    for (auto &dataset : contrast)
      class_counts[dataset.sha1] += dataset.set_size;

  double l = 0;
  const double marginal_motif_prior
      = registration.compute_marginal_motif_prior(present);
  for (auto &contrast : collection)
    for (auto &dataset : contrast) {
      const double class_cond
          = registration.get_class_motif_prior(dataset.sha1, present);
      const double log_class_prior
          = log(registration.get_class_prior(dataset.sha1));

      double motif_count = 0;
#pragma omp parallel for schedule(static) \
    reduction(+ : l, motif_count) if (DO_PARALLEL)
      for (size_t i = 0; i < dataset.set_size; i++) {
        posterior_t res = posterior_atleast_one(dataset.sequences[i], present);
        double p = res.posterior;
        double x = log_class_prior + log(p * class_cond / marginal_motif_prior
                                         + (1 - p) * (1 - class_cond)
                                           / (1 - marginal_motif_prior));
        if (task.measure == Measure::ClassificationLikelihood)
          x += res.log_likelihood;
        if (verbosity >= Verbosity::debug)
          cout << "Sequence " << dataset.sequences[i].definition << " p = " << p
               << " class log likelihood = " << x << " exp -> " << exp(x)
               << endl;
        motif_count += p;
        l += x;
      }
      motif_counts[dataset.sha1] += motif_count;
      if (verbosity >= Verbosity::debug)
        cout << "Data::Set " << dataset.path << " l = " << l << endl;
    }

  // add pseudo counts
  double class_z = 0;
  for (auto &x : class_counts)
    class_z += x.second += 2;
  for (auto &x : motif_counts)
    x.second++;

  // reestimate
  for (auto &x : motif_counts)
    x.second /= class_counts[x.first];
  for (auto &x : class_counts)
    x.second /= class_z;

  double diff = 0;
  score = l;

  // update and compute L1 norm of update
  if (not options.dont_learn_class_prior) {
    if (verbosity >= Verbosity::debug)
      cout << "Updating class prior." << endl;
    for (auto &x : registration.datasets) {
      diff += fabs(x.second.class_prior - class_counts[x.first]);
      x.second.class_prior = class_counts[x.first];
    }
  }
  if (not options.dont_learn_conditional_motif_prior) {
    if (verbosity >= Verbosity::debug)
      cout << "Updating conditional motif priors." << endl;
    for (auto &x : registration.datasets) {
      diff += fabs(x.second.motif_prior[present] - motif_counts[x.first]);
      x.second.motif_prior[present] = motif_counts[x.first];
    }
  }

  if (verbosity >= Verbosity::info) {
    cout << "Re-estimation, L1 norm of class params change  " << diff << endl;
    if (diff < options.termination.gamma_tolerance)
      cout << "L1 norm change criterion                       OK" << endl;
  }

  if (diff < options.termination.gamma_tolerance)
    return true;
  else
    return false;
}

bool HMM::perform_training_iteration_reestimation(
    const Data::Collection &collection, const Training::Task &task,
    const Options::HMM &options, double &score) {
  bool done = false;

  HMM previous(*this);
  double new_score = reestimationIteration(collection, task, options);

  double gamma = norml1(previous.emission - emission)
                 + norml1(previous.transition - transition);
  double delta = -(new_score - score) / new_score;

  if (verbosity >= Verbosity::info) {
    cout << "Re-estimation, L1 norm of parameter change     " << gamma << endl;
    if (gamma < options.termination.gamma_tolerance)
      cout << "L1 norm change criterion                       OK" << endl;

    if (verbosity >= Verbosity::debug
        or gamma < options.termination.gamma_tolerance) {
      cout << "Re-estimation, relative score difference       " << delta
           << endl;
      if (delta < options.termination.delta_tolerance)
        cout << "Relative score criterion                       OK" << endl;
    }
  }

  score = new_score;

  if (gamma < options.termination.gamma_tolerance
      and delta < options.termination.delta_tolerance)
    done = true;

  return done;
}

bool HMM::perform_training_iteration_gradient(
    const Data::Collection &collection, const Training::Task &task,
    const Options::HMM &options, int &center, double &score) {
  if (verbosity >= Verbosity::verbose)
    cerr << "HMM::perform_training_iteration_gradient" << endl;

  double past_score = score;
  bool done = false;

  if (store_intermediate)
    serialize(cerr, options.exec_info);

  if (verbosity >= Verbosity::verbose)
    cerr << "Current parameters: " << *this << endl;

  Timer timer;
  double previous_score;
  Gradient gradient
      = compute_gradient(collection, previous_score, task, options.weighting);
  double gradient_comp_time = timer.tock();
  if (options.timing_information)
    cerr << "Gradient computation time: "
            + time_to_pretty_string(gradient_comp_time) << endl;

  if (verbosity >= Verbosity::info)
    for (size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if (groups[group_idx].name == task.motif_name)
        if (is_motif_group(group_idx))
          cout << "Motif gradient                                 "
               << groups[group_idx].name << ":"
               << get_group_consensus(gradient.emission, group_idx)
               << endl;

  if (verbosity >= Verbosity::verbose)
    cout << "Past score = " << past_score << endl
         << "Previous score = " << previous_score << endl;
  if (verbosity >= Verbosity::verbose)
    cout << "The transition gradient is : " << gradient.transition << endl
         << "The emission gradient is : " << gradient.emission << endl;

  double gradient_norm = sqrt(scalar_product(gradient, gradient));
  if (verbosity >= Verbosity::verbose)
    cout << "Gradient norm = " << gradient_norm << endl;

  double max_abs_gradient_component = 0;
  for (size_t i = 0; i < gradient.emission.size1(); i++)
    for (size_t j = 0; j < gradient.emission.size2(); j++)
      max_abs_gradient_component = max<double>(max_abs_gradient_component,
                                               fabs(gradient.emission(i, j)));
  for (size_t i = 0; i < gradient.transition.size1(); i++)
    for (size_t j = 0; j < gradient.transition.size2(); j++)
      max_abs_gradient_component = max<double>(max_abs_gradient_component,
                                               fabs(gradient.transition(i, j)));

  double new_score = previous_score;
  HMM candidate = *this;

  const double min_score = 0;
  const double min_gradient = 0;
  // const double min_score = 1e-10;
  // const double min_gradient = 1e-15;
  if (previous_score < min_score and gradient_norm < min_gradient) {
    if (verbosity >= Verbosity::info)
      cout << "Skipping line search. Score = " << previous_score
           << " Gradient norm = " << gradient_norm << endl;
  } else if (max_abs_gradient_component == 0) {
    if (verbosity >= Verbosity::info)
      cout << "Skipping line search. Maximal absolute gradient component is "
              "zero." << endl;
  } else if (gradient_norm == 0) {
    if (verbosity >= Verbosity::info)
      cout << "Skipping line search. Gradient norm is zero." << endl;
  } else {
    int info;
    pair<double, HMM> res = line_search_more_thuente(
        collection, gradient, previous_score, info, task, options);
    if (info != 1)
      cout << "Line search exit status: " << info << " - "
           << line_search_status(info) << endl;

    new_score = res.first;
    candidate = res.second;

    if (info != 1 and verbosity >= Verbosity::verbose) {
      cout << "Candidate emission = " << candidate.emission << endl;
      cout << "Candidate transition = " << candidate.transition << endl;
      cout << "Candidate emission gradient = " << gradient.emission << endl;
      cout << "Candidate transition gradient = " << gradient.transition << endl;
    }
  }

  if (verbosity >= Verbosity::info)
    cout << "Score                                          " << new_score
         << endl;

  // TODO reconsider building n-step stopping criteria
  // if(options.training_type == Training::Type::hybrid and
  // options.termination.past == 1)
  past_score = previous_score;

  double score_difference = new_score - past_score;
  double relative_score_difference = score_difference / fabs(new_score);

  matrix_t a = candidate.emission;
  matrix_t b = candidate.transition;
  log_transform(a);
  log_transform(b);
  for (size_t i = 0; i < a.size1(); i++)
    for (size_t j = 0; j < a.size2(); j++)
      if (std::isinf(a(i, j)))
        a(i, j) = 0;
  for (size_t i = 0; i < b.size1(); i++)
    for (size_t j = 0; j < b.size2(); j++)
      if (std::isinf(b(i, j)))
        b(i, j) = 0;

  double xnorm = norml2(a) + norml2(b);
  double gnorm = norml2(gradient.emission) + norml2(gradient.transition);

  if (verbosity >= Verbosity::info)
    cout << "Gradient learning, relative score difference   "
         << relative_score_difference << endl;
  if (verbosity >= Verbosity::verbose) {
    cout << "Score difference = " << score_difference << endl;
    cout << "Gradient norm = " << gnorm << endl;
    cout << "Parameter norm = " << xnorm << endl;
  }

  if (score_difference < 0) {
    cout << "Warning: negative score difference during gradient learning "
            "iteration! Discarding gradient iteration." << endl;
    done = true;
  }

  if (not done) {
    bool rel_score_criterion = relative_score_difference
                               < options.termination.delta_tolerance;
    bool grad_norm_criterion = gnorm < options.termination.epsilon_tolerance
                                       * max<double>(1, xnorm);
    if (rel_score_criterion)
      cout << "Relative score criterion                       OK" << endl;
    if (grad_norm_criterion)
      cout << "Gradient norm criterion                        OK" << endl;
    done = rel_score_criterion or grad_norm_criterion;
    score = new_score;
    *this = candidate;
  } else {
    score = previous_score;
  }

  return done;
}

namespace Exception {
namespace HMM {
namespace Learning {
TrainBgTooLate::TrainBgTooLate()
    : runtime_error(
          "Error: train_background() should not be called after motifs have "
          "been added.") {}
GradientNotImplemented::GradientNotImplemented(
    Measures::Continuous::Measure measure)
    : runtime_error("Calculation of " + measure2string(measure)
                    + " gradient is currently not implemented.") {}
}
}
}
