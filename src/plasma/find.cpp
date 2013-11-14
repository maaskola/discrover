/*
 * =====================================================================================
 *
 *       Filename:  find.cpp
 *
 *    Description:  Routines to discover IUPAC based motifs
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include "score.hpp"
#include "count.hpp"
#include "code.hpp"
#include "randomization_test.hpp"
#include <iostream>
#include <set>
#include <thread>
#include "find.hpp"
#include "mask.hpp"
#include "aux.hpp"
#include "align.hpp"
#include "../mcmc/mcmciupac.hpp"
#include "../timer.hpp"

using namespace std;

namespace Seeding {
  Plasma::Plasma(const Options &opt) : options(opt), collection(options.paths, options.revcomp, options.n_seq), index_ready(false), needs_rebuilding(false) {
    if(options.verbosity >= Verbosity::verbose)
      cerr << "Data loaded - constructor 1." << endl;

    if(options.verbosity >= Verbosity::debug) {
      cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
      cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
      cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
    }


    // TODO reenable this check further down in the code
    // if((options.objective.measure == Measures::Discrete::Measure::signal_frequency or options.objective.measure == Measures::Discrete::Measure::control_frequency)
    // j    and options.plasma.degeneracies.empty() and options.plasma.rel_degeneracy == 1)
    // j  options.plasma.rel_degeneracy = 0.2;

    if(options.plasma.degeneracies.empty() or
        find_if(begin(options.plasma.degeneracies), end(options.plasma.degeneracies), [](size_t a) { return(a!=0); }) != end(options.plasma.degeneracies))
      needs_rebuilding = true;
  }

  Plasma::Plasma(const DataCollection &collection_, const Options &opt) : options(opt), collection(collection_), index_ready(false), needs_rebuilding(false) {
    if(options.verbosity >= Verbosity::verbose)
      cerr << "Data loaded - constructor 2." << endl;

    if(options.verbosity >= Verbosity::debug) {
      cout << "motif_specifications:"; for(auto &x: options.motif_specifications) cout << " " << x; cout << endl;
      cout << "paths:"; for(auto &x: options.paths) cout << " " << x; cout << endl;
      cout << "objectives:"; for(auto &x: options.objectives) cout << " " << x; cout << endl;
    }

    // TODO reenable this check further down in the code
    // if((options.objective.measure == Measures::Discrete::Measure::signal_frequency or options.objective.measure == Measures::Discrete::Measure::control_frequency)
    //     and options.plasma.degeneracies.empty() and options.plasma.rel_degeneracy == 1)
    //   options.plasma.rel_degeneracy = 0.2;

    if(options.plasma.degeneracies.empty() or
        find_if(begin(options.plasma.degeneracies), end(options.plasma.degeneracies), [](size_t a) { return(a!=0); }) != end(options.plasma.degeneracies))
      needs_rebuilding = true;
  }

  void report(ostream &os, const Objective &objective, const string &motif, const DataCollection &collection, const Options &options) {
    Result result(objective);
    result.motif = motif;
    result.counts = count_motif(collection, motif, options);
   // = {motif, 0, 0, counts};
    // Plasma::Result result = {motif, 0, 0, counts};
    report(os, result, collection, options);
  }

  void report(ostream &os, const Result &res, const DataCollection &collection, const Options &options) {
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
    double x = -compute_score(collection, res, options, Measures::Discrete::Measure::CorrectedLogpGtest);
    for(auto &objective: options.objectives)
      os << "Objective: " << Specification::to_string(objective) << "                            " << compute_score(collection, res.counts, options, objective, res.motif.size(), motif_degeneracy(res.motif)) << endl;
    os << "Score                             " << compute_score(collection, res, options) << endl;
    os << "Delta frequency                   " << compute_score(collection, res, options, Measures::Discrete::Measure::DeltaFrequency) << endl;
    os << "Matthew's correlation coefficient " << compute_score(collection, res, options, Measures::Discrete::Measure::MatthewsCorrelationCoefficient) << endl;
    os << "Mutual information [bit]          " << compute_score(collection, res, no_ps, Measures::Discrete::Measure::MutualInformation) << endl;
    os << "Mutual information (pscnt) [bit]  " << compute_score(collection, res, options, Measures::Discrete::Measure::MutualInformation) << endl;
    os << "Expected mutual information [bit] " << compute_score(collection, res, options, Measures::Discrete::Measure::MutualInformation, true) << endl;
    os << "Variance mutual information [bit] " << compute_score(collection, res, options, Measures::Discrete::Measure::VarianceMutualInformation, true) << endl;
    os << "Sd mutual information [bit]       " << sqrt(compute_score(collection, res, options, Measures::Discrete::Measure::VarianceMutualInformation, true)) << endl;
    os << "Z mutual information [bit]        " << compute_score(collection, res, options, Measures::Discrete::Measure::MutualInformation, true) / sqrt(compute_score(collection, res, options, Measures::Discrete::Measure::VarianceMutualInformation, true)) << endl;
    os << "G-test                            " << compute_score(collection, res, options, Measures::Discrete::Measure::Gtest) << endl;
    os << "Uncorrected log-P(G-test)         " << -compute_score(collection, res, options, Measures::Discrete::Measure::LogpGtest) << endl;
    os << "Corrected log-P(G-test)           " << x << (x > 0 ? "      WARNING: greater zero!" : "") << endl;
    size_t i = 0;
    for(auto &series: collection)
      for(auto &set: series) {
        double x = res.counts[i++];
        double z = options.word_stats ? set.seq_size : set.set_size;
        os << "Occurrence statistics             " << set.path << " " << x << " / " << z << " = " << (x/z) << endl;
      }
    if(options.dump_viterbi)
      viterbi_dump(res.motif, collection, cerr, options);
  }

  bool admissible(const std::string &word, const score_map_t &previous, size_t max_degeneracy) {
    return (previous.find(word) == end(previous)) and (motif_degeneracy(word) <= max_degeneracy);
  }

  std::list<std::string> modifications(const std::string &original_word, const std::string &word, size_t position) {
    std::list<std::string> l;
    char c = tolower(original_word[position]);
    std::string replacements;
    switch(c) {
      case 'a':
        replacements = "amrwvhdn";
        break;
      case 'c':
        replacements = "msyvhbn";
        break;
      case 'g':
        replacements = "grskvbdn";
        break;
      case 't':
      case 'u':
        replacements = "twykhbdn";
        break;
      case 'n':
        replacements = "acgtmrwsykvhbdn";
        break;
      default:
        throw("Error: unrecognized letter.");
    }
    for(auto &r: replacements)
      if(r != word[position]) {
        string w(word);
        w[position] = r;
        l.push_back(w);
      }
    return(l);
  }

  Results Plasma::find_seeds(size_t length, const Objective &objective, Algorithm algorithm) {

    size_t max_degeneracy;
    set<size_t> degeneracies;
    for(auto &d: options.plasma.degeneracies)
      degeneracies.insert(d);
    if(degeneracies.empty())
      max_degeneracy = 3 * length;
    else
      max_degeneracy = *degeneracies.rbegin();
    max_degeneracy = min<size_t>(max_degeneracy, 3 * length * options.plasma.rel_degeneracy);
    if(options.plasma.per_degeneracy) {
      for(size_t i = 0; i <= max_degeneracy; i++)
        degeneracies.insert(i);
    }

    if(needs_rebuilding and max_degeneracy > 0)
      rebuild_index();

    Results plasma_results;
    if((algorithm & Algorithm::Plasma) == Algorithm::Plasma)
      plasma_results = find_breadth(length, objective, max_degeneracy, degeneracies);

    Results fire_results;
    if((algorithm & Algorithm::FIRE) == Algorithm::FIRE)
      fire_results = find_fire(length, objective, max_degeneracy, degeneracies);

    Results mcmc_results;
    if((algorithm & Algorithm::MCMC) == Algorithm::MCMC)
      mcmc_results = find_mcmc(length, objective, max_degeneracy);

    Results results;
    set<string> motifs;
    for(auto &m: plasma_results)
      if(motifs.find(m.motif) == end(motifs)) {
        motifs.insert(m.motif);
        results.push_back(m);
      }
    for(auto &m: fire_results)
      if(motifs.find(m.motif) == end(motifs)) {
        motifs.insert(m.motif);
        results.push_back(m);
      }
    for(auto &m: mcmc_results)
      if(motifs.find(m.motif) == end(motifs)) {
        motifs.insert(m.motif);
        results.push_back(m);
      }

    return(results);
  }

  /** This executes MCMC to find discriminative IUPAC motifs.
   */
  Results Plasma::find_mcmc(size_t length, const Objective &objective, size_t max_degeneracy) const {
    srand(time(0));
    MCMC::Evaluator<MCMC::Motif> eval(collection, options, objective);
    MCMC::Generator<MCMC::Motif> gen(options, length, max_degeneracy);
    MCMC::MonteCarlo<MCMC::Motif> mcmc(gen, eval, options.verbosity);
    std::vector<double> temperatures;
    std::vector<MCMC::Motif> init;
    double temperature = options.mcmc.temperature;
    for(size_t i = 0; i < options.mcmc.n_parallel; i++) {
      string word;
      for(size_t j = 0; j < length; j++)
        word += "acgt"[rand() % 4];
      init.push_back(word);
      temperatures.push_back(temperature);
      temperature /= 2;
    }
    auto res = mcmc.parallel_tempering(temperatures, init, options.mcmc.max_iter);

    std::string best_motif = "";
    double best_score = -numeric_limits<double>::infinity();
    for(auto &x: res)
      for(auto &y: x)
        if(y.second > best_score) {
          best_motif = y.first;
          best_score = y.second;
        }

    Results results;

    Stats::OccurrenceCounts best_contrast = count_motif(collection, best_motif, options);
    double log_p = -compute_score(collection, best_contrast, options, objective, length, motif_degeneracy(best_motif), Measures::Discrete::Measure::CorrectedLogpGtest);
    Result result(objective);
    result.motif = best_motif;
    result.score = best_score;
    result.log_p = log_p;
    result.counts = best_contrast;
    results.push_back(result);

    if(options.verbosity >= Verbosity::verbose)
      std::cout << "MCMC found: " << best_motif << " " << best_score << " " << log_p << endl;

    return(results);
  }

  rev_map_t Plasma::determine_initial_candidates(size_t length, const Objective &objective, string &best_motif, size_t &n_candidates, double &max_score, Results &results, const set<size_t> &degeneracies) const {
    const size_t degeneracy = 0;
    rev_map_t candidates;

    best_motif = "";
    max_score = -numeric_limits<double>::infinity();

    Timer my_timer;
    if(options.verbosity >= Verbosity::verbose)
      cerr << "Starting to get word counts." << endl;
    hash_map_t word_counts = get_word_counts(collection, length, options);
    if(options.measure_runtime) {
      cerr << "Got words for length " << length << " in " << my_timer.tock() << " \u00b5s." << endl;
      my_timer.tick();
    }

    //    if(options.verbosity >= Verbosity::debug)
    //      print_counts(word_counts);


    // TODO: instead of just using the N words, the FIRE approach uses randomization test for significance,
    // and stops once 10 consecutive words have failed the test
    for(auto &iter: word_counts) {
      if(options.verbosity >= Verbosity::debug)
        cout << "Candidate " << iter.first << endl;
      double score = compute_score(collection, iter.second, options, objective, length, degeneracy);

      if(options.verbosity >= Verbosity::debug)
        cout << "score = " << score << endl;
      if(score > max_score) {
        max_score = score;
        best_motif = iter.first;
        if(options.verbosity >= Verbosity::debug)
          cout << "motif = " << best_motif << " score = " << score << " " << iter.second << endl;
      }
      if(candidates.empty() or score > candidates.rbegin()->first or n_candidates < options.plasma.max_candidates) {
        candidates.insert({score, iter.first});
        n_candidates++;
        if(options.candidate_selection == CandidateSelection::TopN)
          if(n_candidates > options.plasma.max_candidates) {
            candidates.erase(--end(candidates));
            n_candidates--;
          }
      }
    }

    if(options.candidate_selection == CandidateSelection::RandomizationTest) {
      size_t consecutive_failures = 0;

      rev_map_t candidates_ = candidates;
      rev_map_t candidates = rev_map_t();
      n_candidates = 0;

      for(auto &candidate: candidates_) {
        double score = candidate.first;

        if(options.verbosity >= Verbosity::debug)
          cout << "Candidate " << candidate.second << endl;
        if(options.verbosity >= Verbosity::debug)
          cout << "score = " << score << endl;

        size_t nr_of_tests_to_do = 1 << (2 * length);
        if(options.revcomp)
          nr_of_tests_to_do /= 2;
        bool test_res = randomization_test(collection, word_counts[candidate.second], nr_of_tests_to_do, score, options, objective, length, degeneracy);
        if(test_res)
          consecutive_failures = 0;
        else {
          consecutive_failures++;
          if(consecutive_failures >= options.fire.nr_rand_tests)
            break;
        }

        candidates.insert(candidate);
        n_candidates++;
      }
      // TODO delete the worst options.fire.nr_rand_tests candidates
      if(options.verbosity >= Verbosity::verbose)
        cout << "FIRE got " << n_candidates << " candidates. " << endl;
    }

    if(degeneracies.find(degeneracy) != end(degeneracies)) {
      Stats::OccurrenceCounts best_contrast = count_motif(collection, best_motif, options);
      double log_p = -compute_score(collection, best_contrast, options, objective, length, degeneracy, Measures::Discrete::Measure::CorrectedLogpGtest);
      Result result(objective);
      result.motif = best_motif;
      result.score = max_score;
      result.log_p = log_p;
      result.counts = best_contrast;
      //      Result result = {best_motif, max_score, log_p, best_contrast};
      // cout << "special0: " << best_motif << " " << max_score << " " << best_contrast << endl;
      // report(result, collection, options);
      results.push_back(result);
    }
    /* for(auto &cand: candidates)
       cout << "Candidate " << Motif(cand.second) << " score = " << cand.first << endl;
       */

    if(options.measure_runtime) {
      cerr << "Initial scoring for length " << length << " took " << my_timer.tock() << " \u00b5s." << endl;
      my_timer.tick();
    }

    return(candidates);
  }




  /** This executes the FIRE algorithm to find discriminative IUPAC motifs.
   */
  Results Plasma::find_fire(size_t length, const Objective &objective, size_t max_degeneracy, const set<size_t> &degeneracies) const {
    Results results;
    if(options.verbosity >= Verbosity::verbose)
      cout << "Finding motif of length " << length << " using the FIRE approach with top " << options.plasma.max_candidates << " candidates by " << measure2string(objective.measure) << "." << endl;

    string best_motif;
    size_t n_candidates = 0;
    double max_score;

    rev_map_t candidates = determine_initial_candidates(length, objective, best_motif, n_candidates, max_score, results, degeneracies);

    // add undetermined nucleotides on each side
    for(auto &candidate: candidates) {
      for(size_t i = 0; i < options.fire.nucleotides_5prime; i++)
        candidate.second = "n" + candidate.second;
      for(size_t i = 0; i < options.fire.nucleotides_3prime; i++)
        candidate.second = candidate.second + "n";
    }

    score_map_t examined_words;

    if(max_degeneracy > 0) {
      // wait for the index building to finish
      while(not index_ready) {
        if(options.verbosity >= Verbosity::verbose)
          cerr << "Index still not ready." << endl;
        sleep(1);
      }

      for(auto &candidate: candidates) {
        const string original_word = candidate.second;
        const size_t word_len = original_word.size();

        const size_t n_repetitions = 10;

        if(options.verbosity >= Verbosity::debug)
          std::cout << "FIRE considers motif " << original_word << " " << candidate.first << endl;

        for(size_t idx_repetition = 0; idx_repetition < n_repetitions; idx_repetition++) {
          string word = original_word;
          double previous_score = candidate.first;

          bool tried_all_positions = false;

          if(options.verbosity >= Verbosity::debug)
            std::cout << "FIRE starts repetition " << idx_repetition << " for motif " << original_word << endl;

          while(not tried_all_positions) {
            vector<size_t> remaining_positions;
            for(size_t i = 0; i < word.size(); i++)
              remaining_positions.push_back(i);

            score_map_t tried;

            while(not remaining_positions.empty()) {
              if(options.verbosity >= Verbosity::everything) {
                std::cout << "FIRE remaining positions:";
                for(auto &x: remaining_positions)
                  std::cout << " " << x;
                std::cout << std::endl;
              }

              const size_t position_idx = rand() % remaining_positions.size();
              const size_t position = remaining_positions[position_idx];
              remaining_positions.erase(begin(remaining_positions) + position_idx);

              if(options.verbosity >= Verbosity::debug)
                std::cout << "FIRE tries modifications to position " << position << " of motif " << original_word << endl;

              for(auto &modified_word: modifications(original_word, word, position)) {

                // TODO: FIRE enforces a sufficient improvement of the score
                if(not admissible(modified_word, examined_words, max_degeneracy))
                  continue;

                // find occurrences
                vector<size_t> counts_vec;
                if(options.word_stats)
                  counts_vec = index.word_hits_by_file(modified_word);
                else
                  counts_vec = index.seq_hits_by_file(modified_word, options.revcomp);
                Stats::OccurrenceCounts counts(counts_vec.size());
                for(size_t i = 0; i < counts_vec.size(); i++)
                  counts[i] = counts_vec[i];
                if(options.word_stats and options.revcomp) {
                  counts_vec = index.word_hits_by_file(reverse_complement(modified_word));
                  for(size_t i = 0; i < counts_vec.size(); i++)
                    counts[i] += counts_vec[i];
                }

                double score = compute_score(collection, counts, options, objective, word_len, motif_degeneracy(modified_word));

                tried[modified_word] = score;
                examined_words[modified_word] = score;

                if(options.verbosity >= Verbosity::debug)
                  std::cout << "FIRE: modified word " << modified_word << " -> " << score << endl;
              }

              double achieved_score = previous_score;
              string best = word;
              for(auto &x: tried)
                if(x.second > achieved_score) {
                  achieved_score = x.second;
                  best = x.first;
                }

              if(best != word) {
                word = best;
                previous_score = achieved_score;
                if(options.verbosity >= Verbosity::verbose)
                  std::cout << "FIRE: increased score: " << word << " " << achieved_score << endl;
                break;
              }
            }

            tried_all_positions = remaining_positions.empty();
          }
        }
      }

      for(auto &x: examined_words)
        if(x.second > max_score) {
          max_score = x.second;
          best_motif = x.first;
        }
    }

    Stats::OccurrenceCounts best_contrast = count_motif(collection, best_motif, options);
    double log_p = -compute_score(collection, best_contrast, options, objective, length, motif_degeneracy(best_motif), Measures::Discrete::Measure::CorrectedLogpGtest);
    Result result(objective);
    result.motif = best_motif;
    result.score = max_score;
    result.log_p = log_p;
    result.counts = best_contrast;
    results.push_back(result);

    if(options.verbosity >= Verbosity::verbose)
      std::cout << "FIRE found: " << best_motif << " " << max_score << " " << log_p << endl;

    return(results);
  }

  /** This executes a progressive algorithm to find discriminative IUPAC motifs.
   * It starts with degeneracy 0 and incrementally allows more degeneracy.
   * For each level of degeneracy the top N generalizations of the motifs of the
   * previous level of degeneracy are determined.
   */
  Results Plasma::find_breadth(size_t length, const Objective &objective, size_t max_degeneracy, const set<size_t> &degeneracies) const {
    Results results;
    if(options.verbosity >= Verbosity::verbose)
      cout << "Finding motif of length " << length << " using top " << options.plasma.max_candidates << " breadth search by " << measure2string(objective.measure) << "." << endl;

//    if(options.verbosity >= Verbosity::debug)
//      os << "set signal / control = " << options.set_sizes.signal.size() << " " << options.set_sizes.control.size() << endl;

    size_t degeneracy = 0;

    string best_motif;
    size_t n_candidates = 0;
    double max_score;

    rev_map_t candidates = determine_initial_candidates(length, objective, best_motif, n_candidates, max_score, results, degeneracies);

    bool best_motif_changed = true;

    if(max_degeneracy > 0) {
      while(not index_ready) {
        if(options.verbosity >= Verbosity::verbose)
          cerr << "Index still not ready." << endl;
        sleep(1);
      }
      Timer my_timer;
      while((not candidates.empty()) and degeneracy < max_degeneracy) {
        degeneracy++;
        if(options.verbosity >= Verbosity::verbose)
          cout << "Next round. We have " << candidates.size() << " candidates." << endl;

        score_map_t propositions; // scores have been computed for these motifs
        // alt_score_map_t propositions; // scores have been computed for these motifs
        // set<string> ignored; // motifs whose score is inferior than one of their specifications

        for(auto &candidate: candidates) {
          double candidate_score = candidate.first;
          string candidate_motif(candidate.second);

          if(options.verbosity >= Verbosity::debug)
            cout << "Considering candidate " << candidate_motif << endl;

          for(auto &code: all_generalizations(candidate_motif)) {
            if(options.verbosity >= Verbosity::debug)
              cout << "Considering generalization " << code << endl;
//            string generalization(code);
            if(options.revcomp) {
              string rc = reverse_complement(code);
              if(lexicographical_compare(begin(code), end(code), begin(rc), end(rc)))
                code = rc;
            }
            // bool previously_considered = true;
            // double generalization_score;

            // retrieve or calculate score
            auto iter = propositions.find(code);
            if(iter != end(propositions))
              // it has been considered before, retrieve score
              iter->second = max<double>(iter->second, candidate_score);
            else
              propositions.insert({code, candidate_score});
          }
        }
        vector<score_map_t::const_iterator> work;
        // for(score_map_t::const_iterator &x: propositions)
        for(score_map_t::const_iterator x = propositions.begin(); x != end(propositions); x++)
          work.push_back(x);
        size_t n = work.size();
        vector<double> scores(work.size());
#pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
          string generalization = work[i]->first;
          // double t1;
          // Timer timer;
          // it hasn't been considered before, calculate score
          vector<size_t> counts_vec;
          if(options.word_stats)
            counts_vec = index.word_hits_by_file(generalization);
          else
            counts_vec = index.seq_hits_by_file(generalization, options.revcomp);
          // t1 = timer.tock();
          Stats::OccurrenceCounts counts(counts_vec.size());
          for(size_t i = 0; i < counts_vec.size(); i++)
            counts[i] = counts_vec[i];
          if(options.word_stats and options.revcomp) {
            counts_vec = index.word_hits_by_file(reverse_complement(generalization));
            for(size_t i = 0; i < counts_vec.size(); i++)
              counts[i] += counts_vec[i];
          }
          // Stats::OccurrenceCounts counts_old = count_motif(collection, generalization, options);
          scores[i] = compute_score(collection, counts, options, objective, length, degeneracy);
        }

        candidates = rev_map_t();
        n_candidates = 0;
        for(size_t i = 0; i < n; i++) {
          double candidate_score = work[i]->second;
          double generalization_score = scores[i];
          // heuristic: may be exact for short enough sequences, with length of sequence it becomes more of a heuristic
          if(generalization_score >= candidate_score) {
            string generalization = work[i]->first;
            if(this->options.verbosity >= Verbosity::debug)
              cout << "ax " << generalization << " " << generalization_score << endl;
            if(candidates.empty() or generalization_score > candidates.rbegin()->first or n_candidates < options.plasma.max_candidates) {
              candidates.insert({generalization_score, generalization});
              n_candidates++;
              if(n_candidates > options.plasma.max_candidates) {
                candidates.erase(--end(candidates));
                n_candidates--;
              }
              if(generalization_score > max_score) {
                if(this->options.verbosity >= Verbosity::debug)
                  cout << "New maximum!" << endl;
                max_score = generalization_score;
                best_motif = work[i]->first;
                best_motif_changed = true;
              }
            }
          }
        }

        if(best_motif_changed and degeneracies.find(degeneracy) != end(degeneracies)) {
          Stats::OccurrenceCounts best_contrast = count_motif(collection, best_motif, options);
          double log_p = -compute_score(collection, best_contrast, options, objective, length, degeneracy, Measures::Discrete::Measure::CorrectedLogpGtest);
          Result result(objective);
          result.motif = best_motif;
          result.score = max_score;
          result.log_p = log_p;
          result.counts = best_contrast;
          results.push_back(result);
          best_motif_changed = false;
        }
        if(options.measure_runtime) {
          cerr << "Degeneracy " << degeneracy << " took " << my_timer.tock() << " \u00b5s." << endl;
          my_timer.tick();
        }
      }

      if(best_motif_changed and degeneracies.find(degeneracy) == end(degeneracies)) {
        Stats::OccurrenceCounts best_contrast = count_motif(collection, best_motif, options);
        double log_p = -compute_score(collection, best_contrast, options, objective, length, degeneracy, Measures::Discrete::Measure::CorrectedLogpGtest);
        Result result(objective);
        result.motif = best_motif;
        result.score = max_score;
        result.log_p = log_p;
        result.counts = best_contrast;
        results.push_back(result);
        best_motif_changed = false;
      }
    }
    return(results);
  }

  /** Finds for all included lengths the n_motifs most discriminative motifs.
   * For each length the most discriminative one is found, its occurrences are
   * masked and the procedure is repeated n_motifs times for each length.
   */
  Results Plasma::find_all(const Specification::Motif &motif, const Objective &objective, size_t n_motifs) const {
    Results results;
    for(auto length: motif.lengths) {
      Plasma plasma(*this);
      for(size_t i = 0; i < n_motifs; i++) {
        Results new_results;
        for(auto &result: plasma.find_seeds(length, objective, options.algorithm)) {
          if(options.verbosity >= Verbosity::verbose)
            cout << "Got results: " << result.motif << " " << result.score << endl;
          bool allowed = not options.strict;
          if(not allowed) {
            if(objective.measure == Measures::Discrete::Measure::CorrectedLogpGtest and result.score > 0) // note that logp are stored as their negative in order to be maximized like the other scores
              allowed = true;
            else {
              double log_p = -compute_score(collection, result, options, Measures::Discrete::Measure::CorrectedLogpGtest);
              if(log_p < 0)
                allowed = true;
            }
          }
          if(allowed)
            new_results.push_back(result);
        }

        if(options.verbosity >= Verbosity::verbose)
          for(auto &result: new_results)
            cout << "Got results 2: " << result.motif << " " << result.score << endl;

        if(new_results.empty())
          break;

        sort(begin(new_results), end(new_results), [](const Result &a, const Result &b) { return(a.log_p <= b.log_p); });

        for(auto &result: new_results) {
          results.push_back(result);
          if(i != options.n_motifs - 1)
            plasma.apply_mask(result);
        }
      }
    }
    return(results);
  }

  /** Finds the most discriminative motif over all included lengths.
   * Occurrences are masked and the procedure is repeated n_motifs times.
   */
  Results Plasma::find_multiple(const Specification::Motif &motif, const Objective &objective, size_t n_motifs) const {
    Plasma plasma(*this);
    Results results;
    for(size_t i = 0; i < n_motifs; i++) {
      Results new_results;
      for(auto length: motif.lengths) {
        for(auto &result: plasma.find_seeds(length, objective, options.algorithm)) {
          if(options.verbosity >= Verbosity::verbose)
            cout << "Got results: " << result.motif << " " << result.score << endl;
          bool allowed = not options.strict;
          if(not allowed) {
            if(objective.measure == Measures::Discrete::Measure::CorrectedLogpGtest and result.score > 0) // note that logp are stored as their negative in order to be maximized like the other scores
              allowed = true;
            else {
              double log_p = -compute_score(collection, result, options, Measures::Discrete::Measure::CorrectedLogpGtest);
              if(log_p < 0)
                allowed = true;
            }
          }
          if(allowed)
            new_results.push_back(result);
        }
      }

      if(options.verbosity >= Verbosity::verbose)
        for(auto &result: new_results)
          cout << "Got results 2: " << result.motif << " " << result.score << endl;

      if(new_results.empty())
        return(results);

      sort(begin(new_results), end(new_results), [](const Result &a, const Result &b) { return(a.log_p <= b.log_p); });

      Result result = *begin(new_results);
      results.push_back(result);
      // Mask if necessary
      if(i != options.n_motifs - 1)
        plasma.apply_mask(result);
    }
    return(results);
  };


  Results Plasma::find(const Specification::Motif &motif_spec, const Objectives &objectives, bool doreport) const {
    if(options.verbosity >= Verbosity::debug)
      cout << "motif_spec = " << motif_spec << " objectives = " << objectives << endl;
    for(auto &objective: objectives)
      if(objective.motif_name == motif_spec.name) {

        Results results;
        if(motif_spec.kind == Specification::Motif::Kind::Seed) {
          Result result(objective);
          result.motif = motif_spec.specification;
          result.counts = count_motif(collection, motif_spec.specification, options);
          result.log_p = -compute_score(collection, result.counts, options, objective, result.motif.length(), motif_degeneracy(result.motif), Measures::Discrete::Measure::CorrectedLogpGtest);
          result.score = compute_score(collection, result.counts, options, objective, result.motif.length(), motif_degeneracy(result.motif));
          results.push_back(result);

          return(results);
        }

        if(options.verbosity >= Verbosity::verbose) {
          cout << "IUPAC regular expression motif finding with libPlasma" << endl;
          cout << "objective.measure = '" << objective.measure << "' objective.motif = '" << objective.motif_name << "' series expr = '";
          for(auto &expr: objective.series_expression)
            cout << expr;
          cout << "'" << endl;
        }

        if(options.verbosity >= Verbosity::debug)
          cout << "motif_spec = " << motif_spec << " objective = " << to_string(objective) << endl;

        Timer t;

        if(options.keep_all)
          for(auto &result: find_all(motif_spec, objective, options.n_motifs)) {
            results.push_back(result);
            if(doreport)
              report(cout, result, collection, options);
          }
        else
          for(auto &result: find_multiple(motif_spec, objective, options.n_motifs)) {
            results.push_back(result);
            if(doreport)
              report(cout, result, collection, options);
          }

        double time = t.tock() * 1e-6;
        if(options.measure_runtime)
          cerr << "Processing took " << time << " seconds." << endl;
        return(results);
      }
    cout << "Error: no objective found for motif specification: " << motif_spec.name << ":" << motif_spec.specification << endl;
    exit(-1);
  }

  void Plasma::apply_mask(const string &motif) {
    ::Seeding::apply_mask(collection, motif, options);
    needs_rebuilding = true;
  }

  void Plasma::apply_mask(const Result &result) {
    apply_mask(result.motif);
  }

  void Plasma::apply_mask(const vector<Result> &results) {
    for(auto &result: results)
      apply_mask(result);
  }

  void Plasma::rebuild_index() {
    index_ready = false;
    thread t([&]() {
      Timer my_timer;
      if(options.verbosity >= Verbosity::verbose)
        cerr << "Starting building of index." << endl;
      index = NucleotideIndex<size_t,size_t>(collection, options.verbosity);
      if(options.measure_runtime)
        cerr << "Built index in " << my_timer.tock() << " \u00b5s." << endl;
      index_ready = true;
      });
    t.detach();
  }

  void viterbi_dump(const string &motif, const DataSet &data_set, ostream &out, const Options &options) {
    out << "# " << data_set.path << " details following" << endl;
    for(auto &seq: data_set) {
      size_t n_sites = 0;
      auto match_iter = begin(seq.sequence);
      while((match_iter = search(match_iter, end(seq.sequence), begin(motif), end(motif), iupac_included)) != end(seq.sequence)) {
        n_sites++;
        match_iter++;
      }
      out << ">" << seq.definition << endl
        << "Viterbi #sites = " << n_sites << " Expected #sites = " << n_sites << " P(#sites>=1) = " << ((n_sites > 0) ? 1 : 0) << " Viterbi log-p = nan" << endl
        //      << iter->to_string() << endl
        << seq.sequence << endl;
//      size_t current = 0;
      auto pos_iter = begin(seq.sequence);
      match_iter = begin(seq.sequence);
      while((match_iter = search(match_iter, end(seq.sequence), begin(motif), end(motif), iupac_included)) != end(seq.sequence)) {
//        size_t pos = distance(iter->begin(), pos_iter);
//        auto pos_iter = match_iter;
//        for(size_t j = current; j < pos; j++)
        while(pos_iter++ != match_iter)
          out << "B";
        for(size_t j = 0; j < motif.length(); j++) {
          out << j;
          match_iter++;
//          pos++;
        }
        pos_iter = match_iter;
//        current = pos;
      }
      while(pos_iter++ != end(seq.sequence))
        out << "B";
      out << endl;
    }
  }



  void viterbi_dump(const string &motif, const DataCollection &collection, ostream &out, const Options &options) {
    for(auto &series: collection)
      for(auto &set: series)
        viterbi_dump(motif, set, out, options);
  }

}

