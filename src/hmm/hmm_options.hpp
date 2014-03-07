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
 *       Filename:  hmm_options.hpp
 *
 *    Description:  Data type for options to train and use HMM
 *
 *        Created:  Thu Aug 4 22:12:31 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef HMM_OPTIONS_HPP
#define HMM_OPTIONS_HPP

#include <string>
#include <vector>
#include "../plasma/options.hpp"
#include "basedefs.hpp"
#include "../executioninformation.hpp"
#include "../plasma/specification.hpp"

enum class Compression {
  none,
  gzip,
  bzip2
};

std::string compression2string(Compression compression);
std::string compression2ending(Compression compression);

enum class Relearning {
  None, // Just add, no relearning
  Reestimation, // Adapt transitions and background emissions
  Full // Re-train both previous and newly added motifs
    // TODO: one might want to add the following modes to train only the most recently added motif:
    //   * using conditional MI
    //   * using the motif's objective
};

struct sampling_options {
  bool do_sampling; // whether to perform Gibbs sampling learning
  int min_size;
  int max_size;
  double temperature;
  double anneal_factor;
  size_t n_indels;
  size_t n_shift;
  size_t n_parallel;
};

struct termination_options {
  size_t max_iter, past;
  double gamma_tolerance;
  double delta_tolerance;
  double epsilon_tolerance;
  bool absolute_improvement;
};

struct LineSearchOptions {
  double mu;
  double eta;
  double delta = 0.66;
  size_t max_steps;
};

struct evaluation_options {
  bool occurrence_table;
  bool summary;
  bool viterbi_path;
  bool ric;
};

struct hmm_options {
  std::vector<Specification::Set> paths;
  Specification::Motifs motif_specifications;
  std::vector<std::string> load_paths;
  std::string label;
  std::vector<std::string> seeds;
  Seeding::Options seeding;
  evaluation_options evaluate;
  size_t n_threads;
  size_t bg_order, n_seq;
  double alpha;
  double contingency_pseudo_count, emission_pseudo_count, transition_pseudo_count;
  size_t n_simulations;
  bool long_names;
  bool class_model;
  bool revcomp;
  bool weighting;
  bool accept_multiple;
  Relearning relearning;
  Compression output_compression;
  size_t left_padding, right_padding;
  bool print_posterior;
  bool timing_information;
  size_t cross_validation_iterations;
  double cross_validation_freq;
  bool store_intermediate; // to write out intermediate parameterizations
  size_t wiggle;
  LineSearchOptions line_search;
  unsigned int random_salt; // seed for the random number generator

  bool learn_class_prior;
  bool learn_conditional_motif_prior;
  double class_prior;
  double conditional_motif_prior1, conditional_motif_prior2;

  Training::Method bg_learning;
  Training::Objectives objectives;

  termination_options termination;

  bool limit_logp; // whether to report min(0,corrected logp) or just corrected logp)
  bool use_mi_to_seed;

  sampling_options sampling;

  double lambda;
  std::vector<std::string> emission_matrix_paths;

  Verbosity verbosity;
  ExecutionInformation exec_info;
};

std::istream &operator>>(std::istream &in, Compression &type);
std::istream &operator>>(std::istream &is, Relearning &relearning);

std::ostream &operator<<(std::ostream &os, const Compression &type);
std::ostream &operator<<(std::ostream &os, const Relearning &relearning);
std::ostream &operator<<(std::ostream &os, const Verbosity &verbosity);
std::ostream &operator<<(std::ostream &os, const LineSearchOptions &options);
std::ostream &operator<<(std::ostream &os, const termination_options &options);
std::ostream &operator<<(std::ostream &os, const sampling_options &options);
std::ostream &operator<<(std::ostream &os, const ExecutionInformation &exec_info);
std::ostream &operator<<(std::ostream &os, const hmm_options &options);

#endif

