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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#include <fstream>
#include "../timer.hpp"
#include "../aux.hpp"
#include "hmm.hpp"

#define DO_PARALLEL 1

std::string line_search_status(int status) {
  std::string msg;
  switch(status) {
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
  return(msg);
}


void HMM::train_background(const Data::Collection &data, const hmm_options &options)
{
  HMM bg_hmm(options.verbosity);

  Training::Task task;
  task.measure = Measure::Likelihood;
  task.targets = {bg_hmm.all_range, bg_hmm.emitting_range};

  bg_hmm.reestimation(data, task, options);

  if(options.verbosity >= Verbosity::debug)
    std::cout << "Resulting model = " << bg_hmm << std::endl << std::endl;

  for(size_t i = 0; i < bg_hmm.transition.size1(); i++)
    for(size_t j = 0; j < bg_hmm.transition.size2(); j++)
      transition(i,j) = bg_hmm.transition(i,j);
  for(size_t i = 0; i < bg_hmm.emission.size1(); i++)
    for(size_t j = 0; j < bg_hmm.emission.size2(); j++)
      emission(i,j) = bg_hmm.emission(i,j);

  if(transition.size1() > bg_hmm.transition.size1()) {
    std::cout << "In train_background(); this should not be called after motifs are already present." << std::endl;
    assert(0);
    exit(-1);
  }

  if(options.verbosity >= Verbosity::debug)
    std::cout << *this << std::endl;
}

void HMM::initialize_bg_with_bw(const Data::Collection &collection, const hmm_options &options)
{
  if(options.verbosity >= Verbosity::info)
    std::cout << "Initializing background with Baum-Welch algorithm." << std::endl;

  hmm_options bg_options = options;
  if(options.verbosity == Verbosity::info)
    bg_options.verbosity = Verbosity::error;

  Timer timer;
  train_background(collection, bg_options);
  double time = timer.tock();

  if(options.timing_information)
    std::cerr << "Background learning: " << time << " micro-seconds" << std::endl;
}

double HMM::train(const Data::Collection &training_data, const Training::Tasks &tasks, const hmm_options &options)
{
  if(options.verbosity >= Verbosity::verbose)
    std::cout << "Model to be evaluated = " << *this << std::endl;
  double delta = 0;
  if(tasks.empty()) {
    if(options.verbosity >= Verbosity::info)
      std::cout << "Not performing training because no training tasks were specified." << std::endl;
  } else {
    bool any_found = false;
    for(auto &group: groups)
      for(auto &task: tasks)
        if(task.motif_name == group.name) {
          any_found = true;
        }
    if(not any_found) {
      if(options.verbosity >= Verbosity::info)
        std::cout << "Skipping training because no motifs specified in the tasks have corresponding states in the HMM." << std::endl;
    } else {
      if(options.verbosity >= Verbosity::info)
        std::cout << "Performing training." << std::endl;
      Timer learning_timer;

      if(options.verbosity >= Verbosity::verbose)
        std::cout << "Registering data sets for class based HMMs." << std::endl;

      for(auto &series: training_data)
        for(auto &data_set: series)
          register_dataset(data_set, (1.0*data_set.set_size)/training_data.set_size, options.conditional_motif_prior1, options.conditional_motif_prior2);

      delta = train_inner(training_data, tasks, options);
      if(options.verbosity >= Verbosity::verbose)
        std::cout << std::endl << "The parameters changed by an L1-norm of " << delta << std::endl;

      double time = learning_timer.tock();
      if(options.timing_information)
        std::cerr << "Learning: " << time << " micro-seconds" << std::endl;

      if(options.verbosity >= Verbosity::debug)
        std::cout << "HMM after training:" << std::endl
          << *this << std::endl;

      std::string store_to = options.label + ".hmm";

      if(options.verbosity >= Verbosity::info)
        std::cout << std::endl << "Parameters stored in " << store_to << std::endl;

      std::ofstream os(store_to.c_str());
      serialize(os, options.exec_info);
    }
  }
  return(delta);
}

double HMM::train_inner(const Data::Collection &data, const Training::Tasks &tasks, const hmm_options &options)
{
  if(tasks.empty())
    return(0);
  else {
    HMM previous(*this);

    if(options.sampling.do_sampling) {
      Training::Task task = tasks[0]; // TODO make it work with multiple tasks
      auto sampling_results = mcmc(data, task, options);
      if(options.verbosity >= Verbosity::info)
        std::cout << "Sampling done. Got " << sampling_results.size() << " results." << std::endl;
      double max = -std::numeric_limits<double>::infinity();
      for(auto y : sampling_results)
        for(auto x : y) {
          if(options.verbosity >= Verbosity::verbose)
            std::cout << "Checking if we have a new maximum: " << x.second << std::endl;
          if(x.second > max) {
            *this = x.first;
            max = x.second;
            if(options.verbosity >= Verbosity::info)
              std::cout << "Found new maximum: " << max << std::endl;
          }
        }
    } else
      iterative_training(data, tasks, options);
    double delta = norml1(previous.emission - emission) + norml1(previous.transition - transition);
    return(delta);
  }
}


// double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Series &data, const Training::Targets &targets, const hmm_options &options)
double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Collection &data, const Training::Targets &targets, const hmm_options &options) const
{
  double log_likel = 0;
  for(auto &series: data)
    for(auto &data_set: series)
      if(data_set.motifs.find("control") == data_set.motifs.end())
        log_likel += BaumWelchIteration(T, E, data_set, targets, options);
  if(verbosity >= Verbosity::debug)
    std::cout << "Done BaumWelchIteration(Collection) log_likel = " << log_likel << std::endl;
  return(log_likel);
}

double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Set &data, const Training::Targets &targets, const hmm_options &options) const
{
  double log_likel = 0;
  // TODO find out if there's a performance benefit to the split-up critical regions.
#pragma omp parallel for reduction(+:log_likel) shared(E, T) if(DO_PARALLEL)
  for(size_t j = 0; j < data.sequences.size(); j++)
    log_likel += BaumWelchIteration(T, E, data.sequences[j], targets);
  if(verbosity >= Verbosity::debug)
    std::cout << "Done BaumWelchIteration(Seqs) log_likel = " << log_likel << std::endl;
  return(log_likel);
}

double HMM::BaumWelchIteration_single(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const
{
  size_t L = s.isequence.size();

  vector_t scale;
  matrix_t f = compute_forward_scaled(s, scale);
  matrix_t b = compute_backward_prescaled(s, scale);

  double log_likel = log_likelihood_from_scale(scale);

  if(not targets.transition.empty()) {
    if(not(T.size1() == n_states and T.size2() == n_states))
      T = zero_matrix(n_states, n_states);

    // for all transitions except the one to the start state
    for(size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if(symbol == empty_symbol)
        for(auto k: targets.transition)
          T(k,start_state) += f(i,k) * transition(k,start_state) * b(i+1,start_state);
      else
        // TODO change the order of the loops
        for(auto k: targets.transition) {
          double f_i_k = f(i,k);
          for(auto suc : succ[k])
            T(k,suc) += f_i_k * transition(k,suc) * emission(suc,symbol) * b(i+1,suc);
        }
    }

    // for the transition to the start_state
    for(auto pre : pred[start_state])
      T(pre,start_state) += f(L,pre) * transition(pre,start_state) * b(L+1,start_state);
  }

  if(not targets.emission.empty()) {
    if(not(E.size1() == n_states and E.size2() == n_emissions))
      E = zero_matrix(n_states, n_emissions);

    for(size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if(symbol != empty_symbol)
        for(auto k: targets.emission)
          E(k, symbol) += f(i+1,k) * b(i+1,k) * scale(i+1);
    }
  }
  return(log_likel);
}

double HMM::BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const
{
  size_t L = s.isequence.size();

  vector_t scale;
  matrix_t f = compute_forward_scaled(s, scale);
  matrix_t b = compute_backward_prescaled(s, scale);

  double log_likel = log_likelihood_from_scale(scale);
  if(verbosity >= Verbosity::debug)
    std::cerr << "log_likel = " << log_likel << std::endl;

  if(not targets.transition.empty()) {
    matrix_t t = zero_matrix(n_states, n_states);

    // for all transitions except the one to the start state
    for(size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if(symbol == empty_symbol)
        for(auto k: targets.transition)
          t(k,start_state) += f(i,k) * transition(k,start_state) * b(i+1,start_state);
      else
        // TODO change the order of the loops
        for(auto k: targets.transition) {
          double f_i_k = f(i,k);
          for(auto suc : succ[k])
            t(k,suc) += f_i_k * transition(k,suc) * emission(suc,symbol) * b(i+1,suc);
        }
    }

    // for the transition to the start_state
    for(auto pre : pred[start_state])
      t(pre,start_state) += f(L,pre) * transition(pre,start_state) * b(L+1,start_state);

    if(verbosity >= Verbosity::debug)
      std::cerr << "t = " << t << std::endl;

#pragma omp critical (update_T)
    T += t;
  }

  if(not targets.emission.empty()) {
    matrix_t e = zero_matrix(n_states, n_emissions);

    for(size_t i = 0; i < L; i++) {
      size_t symbol = s.isequence(i);
      if(symbol != empty_symbol)
        for(auto k: targets.emission)
          e(k, symbol) += f(i+1,k) * b(i+1,k) * scale(i+1);
    }

    if(verbosity >= Verbosity::debug)
      std::cerr << "e = " << e << std::endl;

#pragma omp critical (update_E)
    E += e;
  }
  if(verbosity >= Verbosity::debug)
    std::cout << "Done BaumWelchIteration(Seq) log_likel = " << log_likel << std::endl;
  return(log_likel);
}

double HMM::ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Collection &data, const Training::Targets &targets, const hmm_options &options)
{
  double log_likel = 0;
  for(auto &series: data)
    for(auto &data_set : series)
      log_likel += ViterbiIteration(T, E, data_set, targets, options);
  return(log_likel);
}

double HMM::ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Set &data, const Training::Targets &training_targets, const hmm_options &options)
{
  double log_likel = 0;
#pragma omp parallel for reduction(+:log_likel) shared(E, T) if(DO_PARALLEL)
  for(size_t j = 0; j < data.sequences.size(); j++) {
    StatePath path;
    double cur_log_likel = viterbi(data.sequences[j], path);

    size_t L = data.sequences[j].isequence.size();

    matrix_t t = zero_matrix(n_states, n_states);
    if(not training_targets.transition.empty()) {
      t(start_state, path[0]) += 1;
      for(size_t i = 0; i < L-1; i++)
        t(path[i], path[i+1]) += 1;
      t(path[L-1], start_state) += 1;
    }

    matrix_t e = zero_matrix(n_states, n_emissions);
    if(not training_targets.emission.empty())
      for(size_t i = 0; i < L; i++)
        e(path[i], data.sequences[j].isequence[i]) += 1;

#pragma omp critical
    {
      if(not training_targets.emission.empty())
        E += e;
      if(not training_targets.transition.empty())
        T += t;
    }
    log_likel += cur_log_likel;
  }
  return(log_likel);
}

void HMM::reestimation(const Data::Collection &collection, const Training::Task &task, const hmm_options &options)
{
  Training::Tasks tasks;
  tasks.push_back(task);
  iterative_training(collection, tasks, options);
}

double HMM::reestimationIteration(const Data::Collection &data, const Training::Task &task, const hmm_options &options)
{
  double log_likel = -std::numeric_limits<double>::infinity();
  if(store_intermediate)
    serialize(std::cerr, ExecutionInformation());
  matrix_t T, E;
  if(not task.targets.transition.empty())
    T = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty())
    E = zero_matrix(n_states, n_emissions);

  switch(task.measure) {
    case Measure::Likelihood:
      log_likel = BaumWelchIteration(T, E, data, task.targets, options);
      break;
    case Measure::Viterbi:
      log_likel = ViterbiIteration(T, E, data, task.targets, options);
      break;
    default:
      break;
  }

  M_Step(T, E, task.targets, options);

  if(verbosity >= Verbosity::verbose)
    std::cout << "Done HMM::reestimationIteration log_likel = " << log_likel << std::endl;
  return(log_likel);
};

void HMM::M_Step(const matrix_t &T_, const matrix_t &E_, const Training::Targets &targets, const hmm_options &options)
{
  if(verbosity >= Verbosity::verbose) {
    if(not targets.transition.empty())
      std::cerr << "Before normalization: New transitions = " << T_ << std::endl;
    if(not targets.emission.empty())
      std::cerr << "Before normalization: New emissions = " << E_ << std::endl;
  }

  matrix_t T, E;
  // add pseudo count to expected transitions
  if(not targets.transition.empty()) {
    T = T_;
    for(auto t: targets.transition)
      for(size_t j = 0; j < n_states; j++)
        if(transition(t,j) > 0)
          T(t,j) += options.transition_pseudo_count;
  }

  // add pseudo count to expected emissions
  if(not targets.emission.empty()) {
    E = E_;
    for(auto t: targets.emission)
      for(size_t j = 0; j < n_emissions; j++)
        E(t,j) += options.emission_pseudo_count;
  }

  // normalize transition probabilities
  if(not targets.transition.empty())
    normalize_transition(T);

  // normalize emission probabilities of all but the start state
  if(not targets.emission.empty())
    normalize_emission(E);

  for(auto t: targets.transition)
    for(size_t j = 0; j < n_states; j++)
      transition(t,j) = T(t,j);
  for(auto t: targets.emission)
    for(size_t j = 0; j < n_emissions; j++)
      emission(t,j) = E(t,j);

  if(verbosity >= Verbosity::verbose) {
    if(not targets.transition.empty())
      std::cerr << "New transitions = " << transition << std::endl;
    if(not targets.emission.empty())
      std::cerr << "New emissions = " << emission << std::endl;
  }
}

Gradient HMM::compute_gradient(const Data::Collection &data,
    double &score,
    const Training::Task &task,
    bool weighting) const
{
  if(verbosity >= Verbosity::verbose) {
    std::cerr << "HMM::compute_gradient(Data::Collection)" << std::endl
      << "Task = " << task.motif_name << ":";
    for(auto &expr: task.series_expression)
      std::cerr << expr;
    std::cerr << ":" << task.measure << std::endl;
  }
  Gradient gradient;
  if(not task.targets.transition.empty() and gradient.transition.size1() == 0)
    gradient.transition = zero_matrix(n_states, n_states);
  if(not task.targets.emission.empty() and gradient.emission.size1() == 0)
    gradient.emission = zero_matrix(n_states, n_emissions);

  score = 0;
  double W = 0;
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if(is_motif_group(group_idx))
      if(task.motif_name == groups[group_idx].name)
        for(auto &expr: task.series_expression) {
          auto iter = data.find(expr.series);
          if(iter != end(data)) {
            Gradient g;
            double s = compute_gradient(*iter, g, task, group_idx);
            double w = iter->set_size;
            W += w;
            if(not task.targets.transition.empty())
              gradient.transition += g.transition * expr.sign * (weighting ? w : 1.0);
            if(not task.targets.emission.empty())
              gradient.emission += g.emission * expr.sign * (weighting ? w : 1.0);
            if(verbosity >= Verbosity::verbose)
              std::cerr << "series = " << expr.series << " score = " << s << std::endl;
            score += expr.sign * s * (weighting ? w : 1.0);
          }
        }
  if(weighting) {
    score /= W;
    if(not task.targets.transition.empty())
      gradient.transition /= W;
    if(not task.targets.emission.empty())
      gradient.emission /= W;
  }
  if(verbosity >= Verbosity::verbose)
    std::cerr << "HMM::compute_gradient::end score = " << score << std::endl;
  return(gradient);
}

Gradient HMM::compute_gradient(const Data::Series &data,
    double &score,
    const Training::Task &task) const
{
  if(verbosity >= Verbosity::verbose)
    std::cerr << "HMM::compute_gradient(Data::Series, task)" << std::endl;
  // ObjectiveFunction objFunc = objFunc_for_training(task.method);
  Gradient gradient;
  score = 0;
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if(is_motif_group(group_idx))
      score += compute_gradient(data, gradient, task, group_idx);
  if(verbosity >= Verbosity::verbose)
    std::cerr << "HMM::compute_gradient(Data::Series, task)::end" << std::endl;
  return(gradient);
}

double HMM::compute_gradient(const Data::Series &data,
    Gradient &gradient,
    const Training::Task &task,
    size_t group_idx) const
{
  if(verbosity >= Verbosity::verbose)
    std::cerr << "HMM::compute_gradient(Data::Series, task, group_idx)" << std::endl;
  double score = 0;
  switch(task.measure) {
    case Measure::MutualInformation:
      score = mutual_information_gradient(data, task, group_idx, gradient);
      break;
    case Measure::RankInformation:
      score = rank_information_gradient(data, task, group_idx, gradient);
      break;
    case Measure::ChiSquare:
      score = chi_square_gradient(data, task, group_idx, gradient);
      break;
    case Measure::MatthewsCorrelationCoefficient:
      score = matthews_correlation_coefficient_gradient(data, task, group_idx, gradient);
      break;
    case Measure::LogLikelihoodDifference:
      score = log_likelihood_difference_gradient(data, task, group_idx, gradient);
      break;
    case Measure::DeltaFrequency:
      score = site_frequency_difference_gradient(data, task, group_idx, gradient);
      break;
    case Measure::ClassificationPosterior:
    case Measure::ClassificationLikelihood:
      score = class_likelihood_gradient(data, task, group_idx, gradient);
      break;
    default:
      std::cout << "Calculation of " << measure2string(task.measure) << " gradient is currently not implemented." << std::endl;
      exit(-1);
  }

  if(verbosity >= Verbosity::verbose)
    std::cout << "Gradient calculation yielded a score of " << score << "." << std::endl;
  if(verbosity >= Verbosity::verbose)
    std::cerr << "HMM::compute_gradient(Data::Series, task, group_idx)::end" << std::endl;

  return(score);
}

void HMM::iterative_training(const Data::Collection &data,
    const Training::Tasks &tasks,
    const hmm_options &options)
{
  size_t iteration = 0;

  Training::State ts(tasks.size());

  if(verbosity >= Verbosity::info) {
    std::cout << std::endl << "Iteration                                      " << iteration << std::endl;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if(is_motif_group(group_idx))
        std::cout << "Motif consensus                                " << groups[group_idx].name << ":" << get_group_consensus(group_idx) << std::endl;
  }

  while((iteration++ < options.termination.max_iter or options.termination.max_iter == 0 )
      and perform_training_iteration(data, tasks, options, ts))
    if(verbosity >= Verbosity::info) {
      std::cout << std::endl << "Iteration                                      " << iteration << std::endl;
      // std::cout << "Gradient learning, relative score difference   " << relative_score_difference << std::endl;
      for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
        if(is_motif_group(group_idx))
      std::cout << "Motif consensus                                " << groups[group_idx].name << ":" << get_group_consensus(group_idx) << std::endl;
    }

  if(verbosity >= Verbosity::info) {
    std::cout << std::endl << "Finished after " << iteration << " iterations." << std::endl
      << "Final score" << (tasks.size() > 1 ? "s" : "") << " = ";
    for(auto scores: ts.scores)
      std::cout << " " << *scores.rbegin();
    std::cout << std::endl;
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if(is_motif_group(group_idx))
        std::cout << " " << groups[group_idx].name << ":" << get_group_consensus(group_idx);
    std::cout << std::endl;
  }
}

bool HMM::perform_training_iteration(const Data::Collection &data,
    const Training::Tasks &tasks,
    const hmm_options &options,
    Training::State &ts)
{
  bool done = true;

  size_t task_idx = 0;
  for(auto task: tasks) {
    double score = -std::numeric_limits<double>::infinity();
    if(task.measure != Measure::Undefined) {
      if(store_intermediate)
        serialize(std::cerr, ExecutionInformation());

      if(verbosity >= Verbosity::verbose)
        std::cerr << "Current parameters: " << *this << std::endl;

      // collect the states whose parameters are subject to this task
      Training::Range task_targets(task.targets.transition.size() + task.targets.emission.size());
      copy(begin(task.targets.transition), end(task.targets.transition), begin(task_targets));
      copy(begin(task.targets.emission), end(task.targets.emission), begin(task_targets));
      sort(begin(task_targets), end(task_targets));
      task_targets.resize(std::distance(begin(task_targets), unique(begin(task_targets), end(task_targets))));

      // check if the task states are reachable
      if(find(begin(task_targets), end(task_targets), start_state) == end(task_targets)) {
        // the task states do not include the start state
        // check if the task states are reachable from the non-task states
        double z = 0;
        for(size_t i = 0; i < n_states; i++)
          if(find(begin(task_targets), end(task_targets), i) == end(task_targets))
            for(auto j: task_targets)
              z += transition(i,j);
        if(z == 0) {
          if(verbosity >= Verbosity::info)
            std::cout << "None of the task states are reachable from the non-task states; skipping task." << std::endl;
          continue;
        }
      }

      if(options.termination.past <= ts.scores[task_idx].size())
        score = *(ts.scores[task_idx].rbegin() + options.termination.past - 1);

      if(Training::measure2method(task.measure) == Training::Method::Gradient)
        done = perform_training_iteration_gradient(data, task, options, ts.center, score) and done;

      if(Training::measure2method(task.measure) == Training::Method::Reestimation)
        done = perform_training_iteration_reestimation(data, task, options, score) and done;

      if((task.measure == Measure::ClassificationPosterior or task.measure == Measure::ClassificationLikelihood)
          and (options.learn_class_prior or options.learn_conditional_motif_prior))
        done = reestimate_class_parameters(data, task, options, score) and done;

    }
    ts.scores[task_idx++].push_back(score);
  }

  if(store_intermediate)
    serialize(std::cerr, ExecutionInformation());

  return(not done);
}

bool HMM::reestimate_class_parameters(const Data::Collection &collection,
    const Training::Task &task,
    const hmm_options &options,
    double &score)
{
  if(verbosity >= Verbosity::debug)
    std::cout << "HMM::reestimate_class_parameters(Data::Collection, ...)" << std::endl;

  std::unordered_map<std::string, double> class_counts;
  std::unordered_map<std::string, std::map<size_t,double>> motif_counts;
  for(auto &series: collection)
    for(auto &data: series)
      class_counts[data.sha1] += data.set_size;

  double l = 0;
  for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
    if(groups[group_idx].kind == Group::Kind::Motif) {
      const double marginal_motif_prior = compute_marginal_motif_prior(group_idx);
      for(auto &series: collection)
        for(auto &data: series) {
          const double class_cond = get_class_motif_prior(data.sha1, group_idx);
          const double log_class_prior = log(get_class_prior(data.sha1));

          double motif_count = 0;
#pragma omp parallel for schedule(static) reduction(+:l,motif_count) if(DO_PARALLEL)
          for(size_t i = 0; i < data.set_size; i++) {
            // TODO FIX ABSENT
            size_t present_mask = 1 << group_idx;
            size_t absent_mask = 0;
            posterior_t res = posterior_atleast_one(data.sequences[i], present_mask, absent_mask);
            double p = res.posterior;
            double x = log_class_prior + log(p * class_cond / marginal_motif_prior + (1-p) * (1-class_cond) / (1-marginal_motif_prior));
            if(task.measure == Measure::ClassificationLikelihood)
              x += res.log_likelihood;
            if(verbosity >= Verbosity::debug)
              std::cout << "Sequence " << data.sequences[i].definition << " p = " << p << " class log likelihood = " << x << " exp -> " << exp(x) << std::endl;
            motif_count += p;
            l += x;
          }
          motif_counts[data.sha1][group_idx] += motif_count;
          if(verbosity >= Verbosity::debug)
            std::cout << "Data::Set " << data.path << " l = " << l << std::endl;
        }
    }

  // add pseudo counts
  double class_z = 0;
  for(auto &x: class_counts) class_z += x.second += 2;
  for(auto &x: motif_counts) for(auto &y: x.second) y.second++;

  // reestimate
  for(auto &x: motif_counts) for(auto &y: x.second) y.second /= class_counts[x.first];
  for(auto &x: class_counts) x.second /= class_z;

  double diff = 0;
  score = l;

  // update and compute L1 norm of update
  if(options.learn_class_prior) {
    if(verbosity >= Verbosity::debug)
      std::cout << "Updating class prior." << std::endl;
    for(auto &x: registered_datasets) {
      diff += fabs(x.second.class_prior - class_counts[x.first]);
      x.second.class_prior = class_counts[x.first];
    }
  }
  if(options.learn_conditional_motif_prior) {
    if(verbosity >= Verbosity::debug)
      std::cout << "Updating conditional motif priors." << std::endl;
    for(auto &x: registered_datasets)
      for(auto &y: motif_counts[x.first]) {
        diff += fabs(x.second.motif_prior[y.first] - y.second);
        x.second.motif_prior[y.first] = y.second;
      }
  }

  if(verbosity >= Verbosity::info) {
    std::cout << "Re-estimation, L1 norm of class params change  " << diff << std::endl;
    if(diff < options.termination.gamma_tolerance)
      std::cout << "L1 norm change criterion                       OK" << std::endl;
  }

  if(diff < options.termination.gamma_tolerance)
    return(true);
  else
    return(false);
}

bool HMM::perform_training_iteration_reestimation(const Data::Collection &data,
    const Training::Task &task,
    const hmm_options &options,
    double &score)
{
  bool done = false;

  HMM previous(*this);
  double new_score = reestimationIteration(data, task, options);

  double gamma = norml1(previous.emission - emission) + norml1(previous.transition - transition);
  double delta = - (new_score - score) / new_score;

  if(verbosity >= Verbosity::info) {
    std::cout << "Re-estimation, L1 norm of parameter change     " << gamma << std::endl;
    if(gamma < options.termination.gamma_tolerance)
      std::cout << "L1 norm change criterion                       OK" << std::endl;

    if(verbosity >= Verbosity::debug or gamma < options.termination.gamma_tolerance) {
      std::cout << "Re-estimation, relative score difference       " << delta << std::endl;
      if(delta < options.termination.delta_tolerance)
        std::cout << "Relative score criterion                       OK" << std::endl;
    }
  }

  score = new_score;

  if(gamma < options.termination.gamma_tolerance and delta < options.termination.delta_tolerance)
    done = true;

  return(done);
}

std::string consensus(const matrix_t &m, const Training::Range &states, double threshold=0)
{
  const std::string iupac = "-acmgrsvtwyhkdbn";
  std::string consensus = "";
  for(auto &i: states) {
    char present = 0;
    for(size_t j = 0; j < 4; j++)
      if(m(i,j) >= threshold)
        present |= (1 << j);
    consensus += iupac[present];
  }
  return(consensus);
}

bool HMM::perform_training_iteration_gradient(const Data::Collection &data,
    const Training::Task &task,
    const hmm_options &options,
    int &center,
    double &score)
{
  if(verbosity >= Verbosity::verbose)
    std::cerr <<  "HMM::perform_training_iteration_gradient" << std::endl;

  double past_score = score;
  bool done = false;

  if(store_intermediate)
    serialize(std::cerr, ExecutionInformation());

  if(verbosity >= Verbosity::verbose)
    std::cerr << "Current parameters: " << *this << std::endl;

  Timer timer;
  double previous_score;
  Gradient gradient = compute_gradient(data, previous_score, task, options.weighting);
  double gradient_comp_time = timer.tock();
  if(options.timing_information)
    std::cerr << "Gradient computation time: " << gradient_comp_time << " micro-seconds" << std::endl;

  if(verbosity >= Verbosity::info)
    for(size_t group_idx = 0; group_idx < groups.size(); group_idx++)
      if(groups[group_idx].name == task.motif_name)
        if(is_motif_group(group_idx))
          std::cout << "Motif gradient                                 " << groups[group_idx].name << ":" << consensus(gradient.emission, groups[group_idx].states) << std::endl;

  if(verbosity >= Verbosity::verbose)
    std::cout << "Past score = " << past_score << std::endl
      << "Previous score = " << previous_score << std::endl;
  if(verbosity >= Verbosity::verbose)
    std::cout << "The transition gradient is : " << gradient.transition << std::endl
      << "The emission gradient is : " << gradient.emission << std::endl;

  double gradient_norm = sqrt(scalar_product(gradient, gradient));
  if(verbosity >= Verbosity::verbose)
    std::cout << "Gradient norm = " << gradient_norm << std::endl;

  double max_abs_gradient_component = 0;
  for(size_t i = 0; i < gradient.emission.size1(); i++)
    for(size_t j = 0; j < gradient.emission.size2(); j++)
      max_abs_gradient_component = std::max<double>(max_abs_gradient_component, fabs(gradient.emission(i,j)));
  for(size_t i = 0; i < gradient.transition.size1(); i++)
    for(size_t j = 0; j < gradient.transition.size2(); j++)
      max_abs_gradient_component = std::max<double>(max_abs_gradient_component, fabs(gradient.transition(i,j)));

  double new_score = previous_score;
  HMM candidate = *this;

  const double min_score = 0;
  const double min_gradient = 0;
  // const double min_score = 1e-10;
  // const double min_gradient = 1e-15;
  if(previous_score < min_score and gradient_norm < min_gradient) {
    if(verbosity >= Verbosity::info)
      std::cout << "Skipping line search. Score = " << previous_score << " Gradient norm = " << gradient_norm << std::endl;
  } else if(max_abs_gradient_component == 0) {
    if(verbosity >= Verbosity::info)
      std::cout << "Skipping line search. Maximal absolute gradient component is zero." << std::endl;
  } else if(gradient_norm == 0) {
    if(verbosity >= Verbosity::info)
      std::cout << "Skipping line search. Gradient norm is zero." << std::endl;
  } else {
    int info;
    std::pair<double, HMM> res = line_search_more_thuente(data, gradient, previous_score, info, task, options);
    if(info != 1)
      std::cout << "Line search exit status: " << info << " - " << line_search_status(info) << std::endl;

    new_score = res.first;
    candidate = res.second;

    if(info != 1 and verbosity >= Verbosity::verbose) {
      std::cout << "Candidate emission = " << candidate.emission << std::endl;
      std::cout << "Candidate transition = " << candidate.transition << std::endl;
      std::cout << "Candidate emission gradient = " << gradient.emission << std::endl;
      std::cout << "Candidate transition gradient = " << gradient.transition << std::endl;
    }
  }

  if(verbosity >= Verbosity::info)
    std::cout << "Score                                          " << new_score << std::endl;

  // TODO reconsider building n-step stopping criteria
  // if(options.training_type == Training::Type::hybrid and options.termination.past == 1)
  past_score = previous_score;

  double score_difference = new_score - past_score;
  double relative_score_difference = score_difference / fabs(new_score);

  matrix_t a = candidate.emission;
  matrix_t b = candidate.transition;
  log_transform(a);
  log_transform(b);
  for(size_t i = 0; i < a.size1(); i++)
    for(size_t j = 0; j < a.size2(); j++)
      if(isinf(a(i,j)) != 0)
        a(i,j) = 0;
  for(size_t i = 0; i < b.size1(); i++)
    for(size_t j = 0; j < b.size2(); j++)
      if(isinf(b(i,j)) != 0)
        b(i,j) = 0;

  double xnorm = norml2(a) + norml2(b);
  double gnorm = norml2(gradient.emission) + norml2(gradient.transition);

  if(verbosity >= Verbosity::info)
    std::cout << "Gradient learning, relative score difference   " << relative_score_difference << std::endl;
  if(verbosity >= Verbosity::verbose) {
    std::cout << "Score difference = " << score_difference  << std::endl;
    std::cout << "Gradient norm = " << gnorm << std::endl;
    std::cout << "Parameter norm = " << xnorm << std::endl;
  }

  if(score_difference < 0) {
    std::cout << "Warning: negative score difference during gradient learning iteration! Discarding gradient iteration." << std::endl;
    done = true;
  }

  if(not done) {
    bool rel_score_criterion = relative_score_difference < options.termination.delta_tolerance;
    bool grad_norm_criterion = gnorm < options.termination.epsilon_tolerance * std::max<double>(1, xnorm);
    if(rel_score_criterion)
      std::cout << "Relative score criterion                       OK" << std::endl;
    if(grad_norm_criterion)
      std::cout << "Gradient norm criterion                        OK" << std::endl;
    done = rel_score_criterion or grad_norm_criterion;
    score = new_score;
    *this = candidate;
  } else {
    score = previous_score;
  }

  return(done);
}


