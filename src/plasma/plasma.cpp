/*
 * =====================================================================================
 *
 *       Filename:  plasma.cpp
 *
 *    Description:  Routines to discover IUPAC based motifs
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <random>
#include "score.hpp"
#include "count.hpp"
#include "code.hpp"
#include "randomization_test.hpp"
#include <iostream>
#include <fstream>
#include <set>
#include <thread>
#include "plasma.hpp"
#include "mask.hpp"
#include "../aux.hpp"
#include "align.hpp"
#include "../mcmc/mcmciupac.hpp"
#include "../timer.hpp"
#include "dreme/dreme.hpp"

using namespace std;

namespace Seeding {
Plasma::Plasma(const Options &opt)
    : options(opt),
      collection(options.paths, options.revcomp, options.n_seq),
      index_ready(false),
      needs_rebuilding(false) {
  if (options.verbosity >= Verbosity::verbose)
    cerr << "Data loaded - constructor 1." << endl;

  if (options.verbosity >= Verbosity::debug) {
    cout << "motif_specifications:";
    for (auto &x : options.motif_specifications)
      cout << " " << x;
    cout << endl;
    cout << "paths:";
    for (auto &x : options.paths)
      cout << " " << x;
    cout << endl;
    cout << "objectives:";
    for (auto &x : options.objectives)
      cout << " " << x;
    cout << endl;
  }

  // TODO re-enable this check further down in the code
  // if((options.objective.measure == Measures::Discrete::Measure::signal_frequency
  //     or options.objective.measure == Measures::Discrete::Measure::control_frequency)
  //   and options.plasma.degeneracies.empty() and options.plasma.rel_degeneracy == 1)
  //   options.plasma.rel_degeneracy = 0.2;

  if (options.plasma.degeneracies.empty()
      or find_if(begin(options.plasma.degeneracies),
                 end(options.plasma.degeneracies), [](size_t a) {
           return a != 0;
         }) != end(options.plasma.degeneracies))
    needs_rebuilding = true;
}

Plasma::Plasma(const Collection &collection_, const Options &opt)
    : options(opt),
      collection(collection_),
      index_ready(false),
      needs_rebuilding(false) {
  if (options.verbosity >= Verbosity::verbose)
    cerr << "Data loaded - constructor 2." << endl;

  if (options.verbosity >= Verbosity::debug) {
    cout << "motif_specifications:";
    for (auto &x : options.motif_specifications)
      cout << " " << x;
    cout << endl;
    cout << "paths:";
    for (auto &x : options.paths)
      cout << " " << x;
    cout << endl;
    cout << "objectives:";
    for (auto &x : options.objectives)
      cout << " " << x;
    cout << endl;
  }

  // TODO re-enable this check further down in the code
  // if((options.objective.measure == Measures::Discrete::Measure::signal_frequency
  //     or options.objective.measure == Measures::Discrete::Measure::control_frequency)
  //   and options.plasma.degeneracies.empty() and options.plasma.rel_degeneracy == 1)
  //   options.plasma.rel_degeneracy = 0.2;

  if (options.plasma.degeneracies.empty()
      or find_if(begin(options.plasma.degeneracies),
                 end(options.plasma.degeneracies), [](size_t a) {
           return a != 0;
         }) != end(options.plasma.degeneracies))
    needs_rebuilding = true;
}

void report(ostream &os, const Objective &objective, const string &motif,
            const Collection &collection, const Options &options) {
  Result result(objective);
  result.motif = motif;
  result.counts = count_motif(collection, motif, options);
  report(os, result, collection, options);
}

void report(ostream &os, const Result &res, const Collection &collection,
            const Options &options) {
  // TODO make more comprehensive
  // TODO add word counts
  // TODO print PWM of occurrences

  string rc = reverse_complement(res.motif);
  os << "IUPAC representation              " << res.motif << (options.revcomp ? " / " + rc : "") << endl;
  os << "Regular expression                " << iupac2regex(res.motif) << (options.revcomp ? " / " + iupac2regex(rc) : "") << endl;
  os << "Length                            " << res.motif.length() << endl;
  os << "Information content [bit]         " << information_content(res.motif) << endl;
  os << "Information content [bit / pos]   " << information_content(res.motif) / res.motif.length() << endl;
  os << "Degeneracy                        " << motif_degeneracy(res.motif) << endl;

  Options no_ps = options;
  no_ps.pseudo_count = 0;
  double x = -compute_score(collection, res, options,
                            Measures::Discrete::Measure::CorrectedLogpGtest);
  for (auto &objective : options.objectives)
    os << "Objective: " << Specification::to_string(objective)
       << "                            "
       << compute_score(collection, res.counts, options, objective,
                        res.motif.size(), motif_degeneracy(res.motif)) << endl;
  auto score = compute_score(collection, res, options);
  auto dfreq = compute_score(collection, res, options, Measures::Discrete::Measure::DeltaFrequency);
  auto mcc = compute_score(collection, res, options, Measures::Discrete::Measure::MatthewsCorrelationCoefficient);
  auto mi = compute_score(collection, res, no_ps, Measures::Discrete::Measure::MutualInformation);
  auto mi_ps = compute_score(collection, res, options, Measures::Discrete::Measure::MutualInformation);
  auto mi_exp = compute_score(collection, res, options, Measures::Discrete::Measure::MutualInformation, true);
  auto mi_var = compute_score(collection, res, options, Measures::Discrete::Measure::VarianceMutualInformation, true);
  auto mi_sd = sqrt(mi_var);
  auto zscore = mi_exp / mi_sd;
  auto gtest = compute_score(collection, res, options, Measures::Discrete::Measure::Gtest);
  auto logp = -compute_score(collection, res, options, Measures::Discrete::Measure::LogpGtest);
  os << "Score                             " << score << endl;
  os << "Delta frequency                   " << dfreq << endl;
  os << "Matthew's correlation coefficient " << mcc << endl;
  os << "Mutual information [bit]          " << mi << endl;
  os << "Mutual information (pscnt) [bit]  " << mi_ps << endl;
  os << "Expected mutual information [bit] " << mi_exp << endl;
  os << "Variance mutual information [bit] " << mi_var << endl;
  os << "Sd mutual information [bit]       " << mi_sd << endl;
  os << "Z mutual information [bit]        " << zscore << endl;
  os << "G-test                            " << gtest << endl;
  os << "Uncorrected log-P(G-test)         " << logp << endl;
  os << "Corrected log-P(G-test)           " << x
     << (x > 0 ? "      WARNING: greater zero!" : "") << endl;
  size_t i = 0;
  for (auto &contrast : collection)
    for (auto &dataset : contrast) {
      double x = res.counts[i++];
      double z = options.word_stats ? dataset.seq_size : dataset.set_size;
      os << "Occurrence statistics             " << dataset.path << " " << x
         << " / " << z << " = " << (x / z) << endl;
    }
  if (options.dump_viterbi) {
    string viterbi_output = options.label + ".viterbi";
    ofstream ofs(viterbi_output.c_str());
    viterbi_dump(res.motif, collection, ofs, options);
  }
}

Results Plasma::find_seeds(size_t length, const Objective &objective,
                           Algorithm algorithm) {
  size_t max_degeneracy;
  set<size_t> degeneracies;
  for (auto &d : options.plasma.degeneracies)
    degeneracies.insert(d);
  if (degeneracies.empty())
    max_degeneracy = 3 * length;
  else
    max_degeneracy = *degeneracies.rbegin();
  max_degeneracy
      = min<size_t>(max_degeneracy, 3 * length * options.plasma.rel_degeneracy);
  if (options.plasma.per_degeneracy)
    for (size_t i = 0; i <= max_degeneracy; i++)
      degeneracies.insert(i);

  if (needs_rebuilding and max_degeneracy > 0)
    rebuild_index();

  Results plasma_results;
  if ((algorithm & Algorithm::Plasma) == Algorithm::Plasma)
    plasma_results
        = find_plasma(length, objective, max_degeneracy, degeneracies);

  Results external_dreme_results;
  if ((algorithm & Algorithm::ExternalDREME) == Algorithm::ExternalDREME)
    external_dreme_results
        = find_external_dreme(length, objective, max_degeneracy, degeneracies);

  Results mcmc_results;
  if ((algorithm & Algorithm::MCMC) == Algorithm::MCMC)
    mcmc_results = find_mcmc(length, objective, max_degeneracy);

  Results results;
  set<string> motifs;
  for (auto &m : plasma_results)
    if (motifs.find(m.motif) == end(motifs)) {
      motifs.insert(m.motif);
      results.push_back(m);
    }
  for (auto &m : external_dreme_results)
    if (motifs.find(m.motif) == end(motifs)) {
      motifs.insert(m.motif);
      results.push_back(m);
    }
  for (auto &m : mcmc_results)
    if (motifs.find(m.motif) == end(motifs)) {
      motifs.insert(m.motif);
      results.push_back(m);
    }

  return results;
}

/** This executes MCMC to find discriminative IUPAC motifs.
 */
Results Plasma::find_mcmc(size_t length, const Objective &objective,
                          size_t max_degeneracy) const {
  MCMC::Evaluator<MCMC::Motif> eval(collection, options, objective);
  MCMC::Generator<MCMC::Motif> gen(options, length, max_degeneracy);
  MCMC::MonteCarlo<MCMC::Motif> mcmc(gen, eval, options.verbosity);
  std::vector<double> temperatures;
  std::vector<MCMC::Motif> init;
  double temperature = options.mcmc.temperature;
  for (size_t i = 0; i < options.mcmc.n_parallel; i++) {
    init.push_back(gen.generate());
    temperatures.push_back(temperature);
    temperature /= 2;
  }
  auto res = mcmc.parallel_tempering(temperatures, init, options.mcmc.max_iter);

  std::string best_motif = "";
  double best_score = -numeric_limits<double>::infinity();
  for (auto &x : res)
    for (auto &y : x)
      if (y.second > best_score) {
        best_motif = y.first;
        best_score = y.second;
      }

  Results results;

  count_vector_t best_contrast = count_motif(collection, best_motif, options);
  double log_p
      = -compute_score(collection, best_contrast, options, objective, length,
                       motif_degeneracy(best_motif),
                       Measures::Discrete::Measure::CorrectedLogpGtest);
  Result result(objective);
  result.motif = best_motif;
  result.score = best_score;
  result.log_p = log_p;
  result.counts = best_contrast;
  results.push_back(result);

  if (options.verbosity >= Verbosity::verbose)
    std::cout << "MCMC found: " << best_motif << " " << best_score << " "
              << log_p << endl;

  return results;
}

rev_map_t Plasma::determine_initial_candidates(
    size_t length, const Objective &objective, seq_type &best_motif,
    size_t &n_candidates, double &max_score, Results &results,
    const set<size_t> &degeneracies) const {
  const size_t degeneracy = 0;
  rev_map_t candidates;

  best_motif = encode("");
  max_score = -numeric_limits<double>::infinity();

  Timer my_timer;
  if (options.verbosity >= Verbosity::verbose)
    cerr << "Starting to get word counts." << endl;
  hash_map_t word_counts = get_word_counts(collection, length, options);
  if (options.measure_runtime) {
    cerr << "Got words for length " << length << " in " << my_timer.tock()
         << " \u00b5s." << endl;
    my_timer.tick();
  }

  for (auto &iter : word_counts) {
    if (options.verbosity >= Verbosity::debug)
      cout << "Candidate " << decode(iter.first) << endl;
    double score = compute_score(collection, iter.second, options, objective,
                                 length, degeneracy);

    if (options.verbosity >= Verbosity::debug)
      cout << "score = " << score << endl;
    if (score > max_score) {
      max_score = score;
      best_motif = iter.first;
      if (options.verbosity >= Verbosity::debug)
        cout << "motif = " << decode(best_motif) << " score = " << score << " "
             << vec2string(iter.second) << endl;
    }
    if (candidates.empty() or score > candidates.rbegin()->first
        or n_candidates < options.plasma.max_candidates) {
      candidates.insert({score, iter.first});
      n_candidates++;
      if (n_candidates > options.plasma.max_candidates) {
        candidates.erase(--end(candidates));
        n_candidates--;
      }
    }
  }

  if (degeneracies.find(degeneracy) != end(degeneracies)) {
    count_vector_t best_contrast
        = count_motif(collection, decode(best_motif), options);
    double log_p
        = -compute_score(collection, best_contrast, options, objective, length,
                         degeneracy,
                         Measures::Discrete::Measure::CorrectedLogpGtest);
    Result result(objective);
    result.motif = decode(best_motif);
    result.score = max_score;
    result.log_p = log_p;
    result.counts = best_contrast;
    results.push_back(result);
  }

  if (options.measure_runtime) {
    cerr << "Initial scoring for length " << length << " took "
         << my_timer.tock() << " \u00b5s." << endl;
    my_timer.tick();
  }

  return candidates;
}

/** This uses the external program DREME to find discriminative IUPAC motifs.
 */
Results Plasma::find_external_dreme(size_t length, const Objective &objective,
                                    size_t max_degeneracy,
                                    const set<size_t> &degeneracies) const {
  Results results;
  vector<string> paths;
  for (auto &contrast : collection)
    for (auto &dataset : contrast)
      paths.push_back(dataset.path);
  if (paths.size() > 2)
    throw Exception::Dreme::OnlyBinaryContrast();

  string path1 = paths[0];
  string path2 = "";
  if (paths.size() == 2)
    path2 = paths[1];

  auto regexes = Dreme::run(path1, path2, length, length, options.revcomp, 1);

  for (auto &res : regexes) {
    string best_motif = res.first;
    double max_score = res.second;
    count_vector_t best_contrast = count_motif(collection, best_motif, options);
    double log_p
        = -compute_score(collection, best_contrast, options, objective, length,
                         motif_degeneracy(best_motif),
                         Measures::Discrete::Measure::CorrectedLogpGtest);
    Result result(objective);
    result.motif = best_motif;
    result.score = max_score;
    result.log_p = log_p;
    result.counts = best_contrast;
    results.push_back(result);

    if (options.verbosity >= Verbosity::verbose)
      std::cout << "DREME found: " << best_motif << " " << max_score << " "
                << log_p << endl;
  }
  return results;
}

/** This executes a progressive algorithm to find discriminative IUPAC motifs.
 * It starts with degeneracy 0 and incrementally allows more degeneracy.
 * For each level of degeneracy the top N generalizations of the motifs of the
 * previous level of degeneracy are determined.
 */
Results Plasma::find_plasma(size_t length, const Objective &objective,
                            size_t max_degeneracy,
                            const set<size_t> &degeneracies) const {
  Results results;
  if (options.verbosity >= Verbosity::verbose)
    cout << "Finding motif of length " << length << " using top "
         << options.plasma.max_candidates << " breadth search by "
         << measure2string(objective.measure) << "." << endl;

  size_t degeneracy = 0;

  seq_type best_motif;
  size_t n_candidates = 0;
  double max_score;

  rev_map_t candidates = determine_initial_candidates(
      length, objective, best_motif, n_candidates, max_score, results,
      degeneracies);

  bool best_motif_changed = true;

  if (max_degeneracy > 0) {
    while (not index_ready) {
      if (options.verbosity >= Verbosity::verbose)
        cerr << "Index still not ready." << endl;
      sleep(1);
    }
    Timer my_timer;
    while ((not candidates.empty()) and degeneracy < max_degeneracy) {
      degeneracy++;
      if (options.verbosity >= Verbosity::verbose)
        cout << "Next round. We have " << candidates.size() << " candidates."
             << endl;

      score_map_t propositions;  // scores have been computed for these motifs

      for (auto &candidate : candidates) {
        double candidate_score = candidate.first;
        auto candidate_motif(candidate.second);

        if (options.verbosity >= Verbosity::debug)
          cout << "Considering candidate " << decode(candidate_motif) << endl;

        for (auto &code : all_generalizations(candidate_motif)) {
          if (options.verbosity >= Verbosity::debug)
            cout << "Considering generalization " << decode(code) << endl;
          if (options.revcomp) {
            auto rc = iupac_reverse_complement(code);
            if (!lexicographical_compare(begin(code), end(code), begin(rc),
                                         end(rc)))
              code = rc;
          }

          // set or retrieve and update score
          auto iter = propositions.find(code);
          if (iter == end(propositions))
            // it has not been considered before
            // -> set score
            propositions.insert({code, candidate_score});
          else
            // it has been considered before
            // -> retrieve score and possibly update it
            iter->second = max<double>(iter->second, candidate_score);
        }
      }
      vector<score_map_t::const_iterator> work;
      for (score_map_t::const_iterator x = propositions.begin();
           x != end(propositions); x++)
        work.push_back(x);
      const size_t n = work.size();
      vector<double> scores(work.size());
#pragma omp parallel for
      for (size_t i = 0; i < n; i++) {
        auto generalization = work[i]->first;
        count_vector_t counts;
        if (options.word_stats)
          counts = index.word_hits_by_file(generalization, options.revcomp);
        else
          counts = index.seq_hits_by_file(generalization, options.revcomp);
        scores[i] = compute_score(collection, counts, options, objective,
                                  length, degeneracy);
      }

      candidates = rev_map_t();
      n_candidates = 0;
      for (size_t i = 0; i < n; i++) {
        double candidate_score = work[i]->second;
        double generalization_score = scores[i];
        // heuristic: may be exact for short enough sequences, with length of
        // sequence it becomes more of a heuristic
        if (generalization_score >= candidate_score) {
          auto generalization = work[i]->first;
          if (this->options.verbosity >= Verbosity::debug)
            cout << "ax " << decode(generalization) << " "
                 << generalization_score << endl;
          if (candidates.empty()
              or generalization_score > candidates.rbegin()->first
              or n_candidates < options.plasma.max_candidates) {
            candidates.insert({generalization_score, generalization});
            n_candidates++;
            if (n_candidates > options.plasma.max_candidates) {
              candidates.erase(--end(candidates));
              n_candidates--;
            }
            if (generalization_score > max_score) {
              if (this->options.verbosity >= Verbosity::debug)
                cout << "New maximum!" << endl;
              max_score = generalization_score;
              best_motif = work[i]->first;
              best_motif_changed = true;
            }
          }
        }
      }

      if (best_motif_changed
          and degeneracies.find(degeneracy) != end(degeneracies)) {
        count_vector_t best_contrast
            = count_motif(collection, decode(best_motif), options);
        double log_p
            = -compute_score(collection, best_contrast, options, objective,
                             length, degeneracy,
                             Measures::Discrete::Measure::CorrectedLogpGtest);
        Result result(objective);
        result.motif = decode(best_motif);
        result.score = max_score;
        result.log_p = log_p;
        result.counts = best_contrast;
        results.push_back(result);
        best_motif_changed = false;
      }
      if (options.measure_runtime) {
        cerr << "Degeneracy " << degeneracy << " took " << my_timer.tock()
             << " \u00b5s." << endl;
        my_timer.tick();
      }
    }

    if (best_motif_changed
        and degeneracies.find(degeneracy) == end(degeneracies)) {
      count_vector_t best_contrast
          = count_motif(collection, decode(best_motif), options);
      double log_p
          = -compute_score(collection, best_contrast, options, objective,
                           length, degeneracy,
                           Measures::Discrete::Measure::CorrectedLogpGtest);
      Result result(objective);
      result.motif = decode(best_motif);
      result.score = max_score;
      result.log_p = log_p;
      result.counts = best_contrast;
      results.push_back(result);
      best_motif_changed = false;
    }
  }
  return results;
}

/** Finds for all included lengths the n_motifs most discriminative motifs.
 * For each length the most discriminative one is found, its occurrences are
 * masked and the procedure is repeated n_motifs times for each length.
 */
Results Plasma::find_all(const Specification::Motif &motif,
                         const Objective &objective) const {
  Results results;
  for (auto length : motif.lengths) {
    Plasma plasma(*this);
    for (size_t i = 0; i < motif.multiplicity; i++) {
      Results new_results;
      for (auto &result :
           plasma.find_seeds(length, objective, options.algorithm)) {
        if (options.verbosity >= Verbosity::verbose)
          cout << "Got results: " << result.motif << " " << result.score
               << endl;
        bool allowed = not options.strict;
        if (not allowed) {
          if (objective.measure
              == Measures::Discrete::Measure::CorrectedLogpGtest
              and result.score > 0)  // note that logp are stored as their
            // negative in order to be maximized like the
            // other scores
            allowed = true;
          else {
            double log_p = -compute_score(
                               collection, result, options,
                               Measures::Discrete::Measure::CorrectedLogpGtest);
            if (log_p < 0)
              allowed = true;
          }
        }
        if (allowed)
          new_results.push_back(result);
      }

      if (options.verbosity >= Verbosity::verbose)
        for (auto &result : new_results)
          cout << "Got results 2: " << result.motif << " " << result.score
               << endl;

      if (new_results.empty())
        break;

      sort(begin(new_results), end(new_results),
           [](const Result &a, const Result &b) { return a.log_p <= b.log_p; });

      for (auto &result : new_results) {
        results.push_back(result);
        if (i != motif.multiplicity - 1)
          plasma.apply_mask(result);
      }
    }
  }
  return results;
}

/** Finds the most discriminative motif over all included lengths.
 * Occurrences are masked and the procedure is repeated n_motifs times.
 */
Results Plasma::find_multiple(const Specification::Motif &motif,
                              const Objective &objective) const {
  Plasma plasma(*this);
  Results results;
  for (size_t i = 0; i < motif.multiplicity; i++) {
    Results new_results;
    for (auto length : motif.lengths) {
      for (auto &result :
           plasma.find_seeds(length, objective, options.algorithm)) {
        if (options.verbosity >= Verbosity::verbose)
          cout << "Got results: " << result.motif << " " << result.score
               << endl;
        bool allowed = not options.strict;
        if (not allowed) {
          if (objective.measure
              == Measures::Discrete::Measure::CorrectedLogpGtest
              and result.score > 0)
            // note that logp are stored as their negative
            // in order to be maximized like the other scores
            allowed = true;
          else {
            double log_p = -compute_score(
                               collection, result, options,
                               Measures::Discrete::Measure::CorrectedLogpGtest);
            if (log_p < 0)
              allowed = true;
          }
        }
        if (allowed)
          new_results.push_back(result);
      }
    }

    if (options.verbosity >= Verbosity::verbose)
      for (auto &result : new_results)
        cout << "Got results 2: " << result.motif << " " << result.score
             << endl;

    if (new_results.empty())
      return results;

    sort(begin(new_results), end(new_results),
         [](const Result &a, const Result &b) { return a.log_p <= b.log_p; });

    Result result = *begin(new_results);
    results.push_back(result);
    // Mask if necessary
    if (i != motif.multiplicity - 1)
      plasma.apply_mask(result);
  }
  return results;
};

Results Plasma::find_motifs(const Specification::Motif &motif_spec,
                            const Objective &objective, bool doreport) const {
  if (options.verbosity >= Verbosity::debug)
    cout << "motif_spec = " << motif_spec << " objective = " << objective
         << endl;
  if (objective.motif_name != motif_spec.name)
    throw Exception::Plasma::NoObjectiveForMotif(motif_spec.name + ":"
                                                 + motif_spec.specification);

  Results results;
  if (motif_spec.kind == Specification::Motif::Kind::Seed) {
    Result result(objective);
    result.motif = motif_spec.specification;
    result.counts = count_motif(collection, motif_spec.specification, options);
    result.log_p
        = -compute_score(collection, result.counts, options, objective,
                         result.motif.length(), motif_degeneracy(result.motif),
                         Measures::Discrete::Measure::CorrectedLogpGtest);
    result.score
        = compute_score(collection, result.counts, options, objective,
                        result.motif.length(), motif_degeneracy(result.motif));
    results.push_back(result);

    return results;
  }

  if (options.verbosity >= Verbosity::verbose) {
    cout << "IUPAC regular expression motif finding with libPlasma" << endl;
    cout << "objective.measure = '" << objective.measure
         << "' objective.motif = '" << objective.motif_name
         << "' contrast expr = '";
    for (auto &expr : objective.contrast_expression)
      cout << expr;
    cout << "'" << endl;
  }

  if (options.verbosity >= Verbosity::debug)
    cout << "motif_spec = " << motif_spec
         << " objective = " << to_string(objective) << endl;

  Timer t;

  if (not options.only_best)
    for (auto &result : find_all(motif_spec, objective)) {
      results.push_back(result);
      if (doreport)
        report(cout, result, collection, options);
    }
  else
    for (auto &result : find_multiple(motif_spec, objective)) {
      results.push_back(result);
      if (doreport)
        report(cout, result, collection, options);
    }

  double time = t.tock() * 1e-6;
  if (options.measure_runtime)
    cerr << "Processing took " << time << " seconds." << endl;
  return results;
}

void Plasma::apply_mask(const string &motif) {
  ::Seeding::apply_mask(collection, motif, options);
  needs_rebuilding = true;
}

void Plasma::apply_mask(const Result &result) { apply_mask(result.motif); }

void Plasma::apply_mask(const vector<Result> &results) {
  for (auto &result : results)
    apply_mask(result);
}

void Plasma::rebuild_index() {
  index_ready = false;
  thread t([&]() {
    Timer my_timer;
    if (options.verbosity >= Verbosity::verbose)
      cerr << "Starting building of index." << endl;
    index = NucleotideIndex<size_t, size_t>(
        collection, options.allow_iupac_wildcards, options.verbosity);
    if (options.measure_runtime)
      cerr << "Built index in " << my_timer.tock() << " \u00b5s." << endl;
    index_ready = true;
  });
  t.detach();
}

void viterbi_dump(const string &motif, const Set &dataset, ostream &out,
                  const Options &options) {
  out << "# " << dataset.path << " details following" << endl;
  for (auto &seq : dataset) {
    size_t n_sites = 0;
    auto match_iter = begin(seq.sequence);
    while ((match_iter = search(match_iter, end(seq.sequence), begin(motif),
                                end(motif), iupac_included))
           != end(seq.sequence)) {
      n_sites++;
      match_iter++;
    }
    out << ">" << seq.definition << endl << "Viterbi #sites = " << n_sites
        << " Expected #sites = " << n_sites
        << " P(#sites>=1) = " << ((n_sites > 0) ? 1 : 0)
        << " Viterbi log-p = nan" << endl << seq.sequence << endl;
    auto pos_iter = begin(seq.sequence);
    match_iter = begin(seq.sequence);
    while ((match_iter = search(match_iter, end(seq.sequence), begin(motif),
                                end(motif), iupac_included))
           != end(seq.sequence)) {
      while (pos_iter++ != match_iter)
        out << "0";
      for (size_t j = 0; j < motif.length(); j++) {
        out << static_cast<char>('A' + j);
        match_iter++;
      }
      pos_iter = match_iter;
    }
    while (pos_iter++ != end(seq.sequence))
      out << "0";
    out << endl;
  }
}

void viterbi_dump(const string &motif, const Collection &collection,
                  ostream &out, const Options &options) {
  for (auto &contrast : collection)
    for (auto &dataset : contrast)
      viterbi_dump(motif, dataset, out, options);
}

namespace Exception {
namespace Dreme {
const char *OnlyBinaryContrast::what() const noexcept {
  string msg = "Error: DREME can only use binary contrasts to find seeds.";
  return msg.c_str();
}
}
namespace Plasma {
NoObjectiveForMotif::NoObjectiveForMotif(const string &token_)
    : exception(), token(token_) {};
const char *NoObjectiveForMotif::what() const noexcept {
  string msg = "Error: no objective for motif specification: " + token;
  return msg.c_str();
}
}
}
}
