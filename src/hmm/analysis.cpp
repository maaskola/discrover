
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "../aux.hpp"
#include "hmm.hpp"
#include "report.hpp"
#include "../timer.hpp"
#include "../plasma/plasma.hpp"

using namespace std;

// whether to drop models below MICO p-value threshold in multiple motif mode
static const bool drop_below_mico_pvalue_threshold = true;
// MICO p-value threshold to drop models in multiple motif mode
static const double log_pvalue_threshold = log(0.05);
static const Measures::Continuous::Measure filtering_measure
    = Measures::Continuous::Measure::ThresholdedConditionalMutualInformation;

void recreate_symlink(const string &source, const string &target) {
  using namespace boost::filesystem;
  if (exists(target)) {
    if (is_symlink(target))
      remove(target);
    else {
      cout << "Error: file " << target << " exists and is not a symlink. "
           << "File will not be overwritten; aborting." << endl;
      exit(-1);
    }
  }
  create_symlink(absolute(source), target);
}

struct AnalysisResult {
  AnalysisResult(const HMM &model, const Options::HMM &opts)
      : training(),
        training_evaluation(),
        test_evaluation(),
        full_evaluation(),
        model(model),
        options(opts){};
  Training::Result training;
  Evaluator::Result training_evaluation;
  Evaluator::Result test_evaluation;
  Evaluator::Result full_evaluation;
  HMM model;
  Options::HMM options;
  void create_symlinks() const {
    string accepted_label = options.label + ".accepted";
    // create soft-links for his model
    string parameter_path = accepted_label + ".hmm";
    string summary_path = accepted_label + ".summary";
    string table_path = accepted_label + ".table"
                        + compression2ending(options.output_compression);
    string viterbi_path = accepted_label + ".viterbi"
                          + compression2ending(options.output_compression);

    recreate_symlink(training.parameter_file, parameter_path);
    recreate_symlink(full_evaluation.files.summary, summary_path);
    recreate_symlink(full_evaluation.files.table, table_path);
    recreate_symlink(full_evaluation.files.viterbi, viterbi_path);

    vector<string> logo_paths;
    for (auto orig_path : full_evaluation.files.logos) {
      auto here = orig_path.find(options.label);
      if (here != string::npos) {
        auto remainder = orig_path.substr(here + options.label.size());
        string path = accepted_label + remainder;
        recreate_symlink(orig_path, path);
        logo_paths.push_back(path);
      }
    }

    if (options.verbosity >= Verbosity::info) {
      cout << "The results of the accepted model can be found in" << endl;
      cout << parameter_path << endl;
      cout << summary_path << endl;
      cout << table_path << endl;
      cout << viterbi_path << endl;
      for (auto path : logo_paths)
        cout << path << endl;
    }
  }
};

size_t compute_degrees_of_freedom(const Data::Collection &data,
                                  const Options::HMM &options,
                                  Specification::Motif &motif_spec) {
  double df = 0;
  for (auto &contrast : data) {
    bool found = false;
    for (auto &set_spec : options.paths)
      if (set_spec.contrast == contrast.name
          and set_spec.motifs.find(motif_spec.name) != end(set_spec.motifs))
        found = true;
    if (found)
      df += (contrast.sets.size() - 1);
  }
  return df;
}

void check_data(const Data::Collection &collection,
                const Options::HMM &options) {
  if (options.verbosity >= Verbosity::info)
    for (auto &contrast : collection) {
      cout << "Data collection has " << contrast.set_size
           << " sequences with a total size of " << contrast.seq_size
           << " nucleotides." << endl;
      for (auto &dataset : contrast) {
        cout << dataset.path << " has " << dataset.set_size
             << " sequences with a total size of " << dataset.seq_size
             << " nucleotides." << endl;
        if (options.verbosity >= Verbosity::verbose) {
          cout << "The SHA1 of this set is " << dataset.sha1 << endl;
          cout << "The first 3 sequences are:" << endl;
          size_t i = 0;
          for (auto &seq : dataset) {
            cout << ">" << seq.definition << endl << seq.sequence << endl;
            if (++i >= 3)
              break;
          }
        }
        if (options.verbosity >= Verbosity::debug)
          for (auto &seq : dataset)
            cout << ">" << seq.definition << endl << seq.sequence << endl;
      }
    }

  // TODO re-enable warning
  //  if(options.training_method != Training::Method::none and not
  //  is_generative(options.objective.measure) and options.objective.measure !=
  //  Measure::rank_information) {
  //    for(auto &contrast: collection) {
  //      if(contrast.sets.size() < 2) {
  //        cout << "Please note that for discriminative training you need to
  //        specify multiple sequence sets using -f." << endl;
  //        exit(-1);
  //      }
  //    }
  //  }
}

AnalysisResult train_evaluate(HMM &hmm, const Data::Collection &all_data,
                              const Data::Collection &training_data,
                              const Data::Collection &test_data,
                              const Options::HMM &options, bool do_training,
                              bool relearning_phase = false) {
  AnalysisResult result(hmm, options);
  // Define the learning and evaluation tasks
  Training::Tasks eval_tasks = hmm.define_training_tasks(options);

  if (do_training) {
    // define the training tasks
    Training::Tasks learn_tasks = hmm.define_training_tasks(options);

    if (relearning_phase) {
      switch (options.multi_motif.relearning) {
        case Options::MultiMotif::Relearning::None:
          learn_tasks = Training::Tasks();
          break;
        case Options::MultiMotif::Relearning::Full:
          break;
        // case Relearning::Added:
        // TODO implement
        //   break;
        case Options::MultiMotif::Relearning::Reestimation: {
          set<string> contrast_names;
          for (auto &task : eval_tasks)
            for (auto &expr : task)
              contrast_names.insert(expr.contrast);

          learn_tasks = Training::Tasks();
          Training::Task task;
          task.motif_name = "Context";
          task.measure = Measure::Likelihood;
          for (auto &contrast_name : contrast_names)
            task.contrast_expression.push_back({+1, contrast_name});

          task.targets.emission.push_back(1);
          for (size_t i = 0; i < hmm.get_nstates(); i++)
            task.targets.transition.push_back(i);

          if (options.verbosity >= Verbosity::verbose) {
            cout << "Generated generative training targets for re-learning."
                 << endl << "Emissions:";
            for (auto e : task.targets.emission)
              cout << " " << e;
            cout << endl;
            cout << "Transitions:";
            for (auto t : task.targets.transition)
              cout << " " << t;
            cout << endl;
          }
          learn_tasks.push_back(task);
        }
      }
    }

    if (not learn_tasks.empty())
      result.training = hmm.train(training_data, learn_tasks, options);
  }

  Evaluator evaluator(hmm);
  if (test_data.set_size != 0) {
    result.training_evaluation
        = evaluator.report(training_data, "training", eval_tasks, options);
    result.test_evaluation
        = evaluator.report(test_data, "Test", eval_tasks, options);
  }
  result.full_evaluation = evaluator.report(all_data, "", eval_tasks, options);
  result.model = hmm;
  return result;
}

double calc_expected_seq_size(const Data::Collection &collection) {
  return 1.0 * collection.seq_size / collection.set_size;
}

string padding(size_t n, char c = 'n') {
  string s = "";
  for (size_t i = 0; i < n; i++)
    s += c;
  return s;
}

vector<string> generate_wiggle_variants(const string &s, size_t n,
                                        Verbosity verbosity) {
  if (verbosity >= Verbosity::debug)
    cout << "generate_wiggle_variants(" << s << ", " << n << ")" << endl;
  size_t l = s.size();
  vector<string> variants;
  variants.push_back(s);
  for (size_t i = 1; i <= min<size_t>(l, n); i++) {
    variants.push_back(padding(i) + s.substr(0, l - i));
    variants.push_back(s.substr(i, l - i) + padding(i));
  }
  if (verbosity >= Verbosity::debug)
    for (auto &x : variants)
      cout << "Wiggle variant: " << x << endl;
  return variants;
}

void prepare_scoring_model(HMM &model, vector<size_t> &present_groups,
                           vector<size_t> &previous_groups,
                           const Options::HMM &options) {
  if (options.revcomp) {
    if (options.verbosity >= Verbosity::info)
      cout << "Adding reverse complementary motifs." << endl;

    pair<HMM, map<size_t, size_t>> rc = model.add_revcomp_motifs();
    model = rc.first;
    map<size_t, size_t> rc_assoc = rc.second;

    // add reverse-complementary motifs of the present motifs
    auto scoring_present_groups = present_groups;
    for (auto x : present_groups)
      scoring_present_groups.push_back(rc_assoc[x]);
    present_groups = scoring_present_groups;

    // add reverse-complementary motifs of the previous motifs
    auto scoring_previous_groups = previous_groups;
    for (auto x : previous_groups)
      scoring_previous_groups.push_back(rc_assoc[x]);
    previous_groups = scoring_previous_groups;
  }

  if (options.verbosity >= Verbosity::info)
    for (size_t i = 0; i < model.get_ngroups(); i++)
      if (model.is_motif_group(i))
        cout << (find(begin(present_groups), end(present_groups), i)
                     != end(present_groups)
                     ? "+"
                     : "-") << " " << model.get_group_consensus(i) << endl;
}

bool check_pairwise_with_previous_motifs(
    const HMM &model, const Data::Collection &data,
    const vector<size_t> &present_groups, const vector<size_t> &previous_groups,
    const Options::HMM &options, size_t n, size_t df, size_t motif_len,
    Measures::Continuous::Measure measure) {
  if (previous_groups.size() > 1) {
    if (options.verbosity >= Verbosity::info)
      cout << "Filtering motif for redundancy with all individual previous "
              "motifs:" << endl;
    for (auto &previous_group : previous_groups) {
      bool ok = true;

      // TODO: reduce model to only these two motifs?
      auto scoring_model = model;
      auto scoring_present_groups = present_groups;
      vector<size_t> scoring_previous_groups;
      scoring_previous_groups.push_back(previous_group);
      prepare_scoring_model(scoring_model, scoring_present_groups,
                            scoring_previous_groups, options);
      double score = scoring_model.compute_score(data, measure, options,
                                                 scoring_present_groups,
                                                 scoring_previous_groups);
      switch (measure) {
        case Measures::Continuous::Measure::
            ConditionalPairMutualInformationRatio:
          ok = score >= options.multi_motif.residual_ratio;
          break;
        case Measures::Continuous::Measure::ConditionalMutualInformation:
        case Measures::Continuous::Measure::
            ThresholdedConditionalMutualInformation: {
          if (options.verbosity >= Verbosity::info)
            cout << "conditional MICO = " << score << endl;
          double log_pvalue
              = corrected_pvalue(score, n, df, motif_len, Verbosity::verbose);
          if (options.verbosity >= Verbosity::info)
            cout << "log p-value of conditional MICO = " << log_pvalue << endl;
          ok = log_pvalue <= log_pvalue_threshold;
        } break;
        default:
          cout << "Error: measure " << measure
               << " is not implemented for filtering in multiple motif mode."
               << endl;
          exit(-1);
      }

      if (not ok) {
        if (options.verbosity >= Verbosity::info)
          cout << "Motif "
               << model.get_group_consensus(*present_groups.rbegin())
               << " is redundant with previous motif "
               << model.get_group_consensus(previous_group) << endl;
        return false;
      }
    }
  }
  return true;
}

HMM doit(const Data::Collection &all_data,
         const Data::Collection &training_data,
         const Data::Collection &test_data, const Options::HMM &options_) {
  vector<AnalysisResult> results;
  // potentially: one might use the MICO p-value for selection,
  // regardless of the objective function chosen for training!

  Options::HMM options = options_;

  if (options.verbosity >= Verbosity::debug)
    cout << "About to construct HMM." << endl;

  // initialize HMM
  HMM hmm(options.verbosity, options.contingency_pseudo_count);

  // load HMM from file if specified
  size_t n_loaded = 0;
  bool training_necessary = false;
  for (auto &load_path : options.load_paths)
    if (n_loaded++ == 0)
      hmm = HMM(load_path, options.verbosity, options.contingency_pseudo_count);
    else {
      hmm.add_motifs(
          HMM(load_path, options.verbosity, options.contingency_pseudo_count));
      training_necessary = true;
    }

  if (options.verbosity >= Verbosity::verbose)
    cout << "Model after initialization = " << hmm << endl;

  hmm.switch_intermediate(options.store_intermediate);

  // train background
  if (n_loaded == 0 and not options.objectives.empty()) {
    hmm.initialize_bg_with_bw(training_data, options);

    if (options.verbosity >= Verbosity::verbose)
      cout << "Model after background learning = " << hmm << endl;
  }

  double expected_seq_size = calc_expected_seq_size(all_data);
  if (options.verbosity >= Verbosity::info)
    cout << "The average sequence size is " << expected_seq_size << "nt."
         << endl;

  bool run_plasma = false;

  for (auto &spec : options.motif_specifications)
    switch (spec.kind) {
      case Specification::Motif::Kind::File:
        // load emission matrix and add it.
        hmm.add_motif(read_emission(spec.specification), expected_seq_size,
                      options.lambda, spec.name, spec.insertions,
                      options.left_padding, options.right_padding);
        break;
      case Specification::Motif::Kind::Seed:
        // use IUPAC regular expression
        hmm.add_motif(spec.specification, options.alpha, expected_seq_size,
                      options.lambda, spec.name, spec.insertions,
                      options.left_padding, options.right_padding);
        training_necessary = true;
        break;
      case Specification::Motif::Kind::Plasma:
        options.seeding.motif_specifications.push_back(spec);
        training_necessary = true;
        run_plasma = true;
        break;
    }

  if (not run_plasma) {
    if (options.verbosity >= Verbosity::info)
      cout << "No automatic seeds are used." << endl;
    auto result = train_evaluate(hmm, all_data, training_data, test_data,
                                 options, training_necessary);
    results.push_back(result);
  } else {
    if (options.verbosity >= Verbosity::info)
      cout << "Determining seeds automatically." << endl;

    Seeding::Collection collection(training_data);

    if (options.verbosity >= Verbosity::debug) {
      for (auto &contrast : training_data)
        for (auto &dataset : contrast) {
          cerr << "HMM::Contrast " << contrast.name << " set -> motifs:";
          for (auto &m : dataset.motifs)
            cerr << " " << m;
          cerr << endl;
        }
    }
    if (options.verbosity >= Verbosity::debug) {
      for (auto &contrast : collection)
        for (auto &dataset : contrast) {
          cerr << "Seeding::Contrast " << contrast.name << " set -> motifs:";
          for (auto &m : dataset.motifs)
            cerr << " " << m;
          cerr << endl;
        }
    }

    Seeding::Plasma plasma(collection, options.seeding);

    if (options.verbosity >= Verbosity::debug)
      cout << "motif_specs.size() = "
           << plasma.options.motif_specifications.size() << endl;

    // determine matching Plasma objectives
    plasma.options.objectives = Training::corresponding_objectives(
        options.objectives, options.use_mi_to_seed);

    // while there are motif specifications left
    for (auto motif_spec : options.motif_specifications) {
      if (hmm.get_nmotifs() > 0) {
        cout << "Masking plasma data collection." << endl;
        plasma.collection.mask(hmm.compute_mask(training_data));
      }

      auto plasma_objective
          = Seeding::objective_for_motif(plasma.options.objectives, motif_spec);

      if (options.verbosity >= Verbosity::info) {
        cout << "Determining seed " << motif_spec << " with objective function "
             << plasma_objective << endl;
      }

      // find the next set of seeds with Plasma
      Seeding::Results plasma_results
          = plasma.find_motifs(motif_spec, plasma_objective);

      // if no seeds are found, don't try to find more seeds
      if (plasma_results.empty())
        break;
      if (options.verbosity >= Verbosity::debug)
        cout << motif_spec.name << " result.size() = " << plasma_results.size()
             << endl;

      vector<pair<string, HMM>> learned_models;

      // seed and learn HMM parameters independently for each Plasma motif
      for (size_t seed_idx = 0; seed_idx < plasma_results.size(); seed_idx++) {
        string motif = plasma_results[seed_idx].motif;

        cout << "Considering candidate motif " << motif_spec.name << ":"
             << motif << " and training to determine the HMM score." << endl;

        plasma_results[seed_idx].score = -numeric_limits<double>::infinity();

        // consider all wiggle variants
        for (auto variant : generate_wiggle_variants(motif, options.wiggle,
                                                     options.verbosity)) {
          HMM model(hmm);
          if (options.verbosity >= Verbosity::info and options.wiggle > 0)
            cout << "Considering wiggle variant " << variant
                 << " of candidate motif " << motif_spec.name << ":" << motif
                 << " and training to determine the HMM score." << endl;

          if (options.extend > 0) {
            if (options.verbosity >= Verbosity::info)
              cout << "Extending seed by " << options.extend
                   << " nucleotides of N." << endl;
            variant = padding(options.extend) + variant
                      + padding(options.extend);
          }

          model.add_motif(variant, options.alpha, expected_seq_size,
                          options.lambda, motif_spec.name,
                          motif_spec.insertions, options.left_padding,
                          options.right_padding);

          Options::HMM options_(options);
          if (options_.long_names)
            options_.label += "." + variant;
          auto result = train_evaluate(model, all_data, training_data,
                                       test_data, options_, training_necessary);
          results.push_back(result);

          learned_models.push_back(make_pair(variant, model));
        }
      }

      if (not options.multi_motif.accept_multiple) {
        if (learned_models.size() == 1)
          hmm = learned_models.begin()->second;
        else {
          if (options.verbosity >= Verbosity::info)
            cout << "Evaluating learned models.";
          auto best_model = hmm;
          string best_seed = "";
          double best_score = -numeric_limits<double>::infinity();
          for (auto &learned : learned_models) {
            string seed = learned.first;
            auto model = learned.second;
            Training::Tasks tasks = model.define_training_tasks(options);
            // TODO CONSIDER: using compute_score instead of
            // compute_score_all_motifs
            double score = model.compute_score_all_motifs(
                training_data, tasks.begin()->measure, options);
            if (score > best_score) {
              best_score = score;
              best_seed = seed;
              best_model = model;
            }
          }
          if (best_score >= -numeric_limits<double>::infinity()) {
            if (options.verbosity >= Verbosity::info)
              cout << "Accepting seed " << best_seed << " with score "
                   << best_score << endl;
            hmm = best_model;
          } else if (options.verbosity >= Verbosity::info)
            cout << "We did not find an acceptable model." << endl;
        }
      } else {
        bool try_finding_another_motif = true;

        double best_score_for_this_motif = -numeric_limits<double>::infinity();

        vector<size_t> previous_groups;
        while (try_finding_another_motif and not learned_models.empty()) {
          auto data = training_data;

          size_t best_index = 0;
          auto best_model = hmm;
          string best_seed = "";
          bool updated = false;
          double best_log_pvalue = 0;

          size_t index = 0;

          vector<string> below_threshold;

          for (auto &learned : learned_models) {
            string seed = learned.first;
            auto learned_model = learned.second;

            auto model = hmm;
            model.add_motifs(learned_model, false);

            // TODO CONSIDER: re-learn at this point?

            // motif_len, n, and df are needed for MICO p-value computation
            const size_t motif_len = seed.size();
            const double n = data.set_size;  // TODO FIX this with regards to
                                             // pseudo counts and exact
                                             // reference to the relevant
                                             // contrast
            size_t df = compute_degrees_of_freedom(data, options, motif_spec);

            vector<size_t> present_groups
                = {model.get_ngroups()
                   - 1};  // only add the most recently added group

            double mi = 0;
            double pairwise_ok = true;

            if (previous_groups.empty()) {
              mi = model.compute_score(
                  data, Measures::Continuous::Measure::MutualInformation,
                  options, present_groups);
              if (options.verbosity >= Verbosity::info)
                cout << "mi = " << mi << endl;
            } else {
              pairwise_ok = check_pairwise_with_previous_motifs(
                  model, data, present_groups, previous_groups, options, n, df,
                  motif_len, filtering_measure);
              if (pairwise_ok) {
                if (options.verbosity >= Verbosity::info)
                  cout << "Scoring residual information against the "
                          "disjunction of all previous motifs:" << endl;

                auto scoring_model = model;
                auto scoring_present_groups = present_groups;
                auto scoring_previous_groups = previous_groups;
                prepare_scoring_model(scoring_model, scoring_present_groups,
                                      scoring_previous_groups, options);
                switch (filtering_measure) {
                  case Measures::Continuous::Measure::
                      ConditionalMutualInformation:
                  case Measures::Continuous::Measure::
                      ThresholdedConditionalMutualInformation:
                    mi = scoring_model.compute_score(
                        data, filtering_measure, options,
                        scoring_present_groups, scoring_previous_groups);
                    break;
                  default:
                    cout << "Error: measure " << filtering_measure
                         << " is not implemented for filtering in multiple "
                            "motif mode." << endl;
                    exit(-1);
                }

                if (options.verbosity >= Verbosity::info)
                  cout << "residual mutual information = " << mi << endl;
              }
            }

            if (pairwise_ok) {
              double log_pvalue
                  = corrected_pvalue(mi, n, df, motif_len, Verbosity::verbose);

              if (options.verbosity >= Verbosity::info)
                cout << "Score of the model augmented by " << seed
                     << " has a log p-value of " << log_pvalue << endl << endl;

              if ((not drop_below_mico_pvalue_threshold)
                  or log_pvalue <= log_pvalue_threshold) {
                if (log_pvalue < best_log_pvalue) {
                  cout << "This is the currently best model!" << endl << endl;
                  try_finding_another_motif = true;
                  updated = true;
                  best_model = model;
                  best_seed = seed;
                  best_log_pvalue = log_pvalue;
                  best_index = index;
                }
              } else
                below_threshold.push_back(seed);
            } else if (drop_below_mico_pvalue_threshold)
              below_threshold.push_back(seed);
            index++;
          }

          if (not updated) {
            cout << "We did not find an improved model." << endl;
            try_finding_another_motif = false;
          } else {
            hmm = best_model;

            if (options.verbosity >= Verbosity::info)
              cout << "Accepting seed " << best_seed << " with log p-value "
                   << best_log_pvalue << endl;

            Options::HMM options_(options);
            if (options.long_names)
              options.label += "." + best_seed;

            auto result
                = train_evaluate(hmm, all_data, training_data, test_data,
                                 options, training_necessary, true);
            results.push_back(result);

            previous_groups.push_back(hmm.get_ngroups() - 1);
            learned_models.erase(begin(learned_models) + best_index);

            size_t task_idx = 0;  // TODO FIX choose the right task; for now
                                  // assume it is the first one
            double current_score
                = *result.training.state.scores[task_idx].rbegin();
            if (current_score > best_score_for_this_motif) {
              if (options.verbosity >= Verbosity::info)
                cout << "This model improves the best score to "
                     << current_score << endl;
              best_score_for_this_motif = current_score;
              try_finding_another_motif = true;
            } else {
              if (options.verbosity >= Verbosity::info)
                cout << "This model does not improve the best score. Stopping "
                        "training." << endl;
              try_finding_another_motif = false;
            }
          }

          if (options.verbosity >= Verbosity::info) {
            cout << "To be removed:";
            for (auto &mod : below_threshold)
              cout << " " << mod;
            cout << endl;
          }

          auto is_below_threshold =
              [&below_threshold](const pair<string, HMM> &x) {
            return find(begin(below_threshold), end(below_threshold), x.first)
                   != end(below_threshold);
          };

          if (options.verbosity >= Verbosity::info) {
            cout << "Before removal:";
            for (auto &mod : learned_models)
              cout << " " << mod.first;
            cout << endl;
            cout << "Before removal " << learned_models.size() << " remaining."
                 << endl;
          }

          if (drop_below_mico_pvalue_threshold)
            learned_models.erase(
                remove_if(begin(learned_models), end(learned_models),
                          is_below_threshold),
                end(learned_models));

          if (options.verbosity >= Verbosity::info) {
            cout << "After removal:";
            for (auto &mod : learned_models)
              cout << " " << mod.first;
            cout << endl;
            cout << "After removal " << learned_models.size() << " remaining."
                 << endl;
          }
        }

        // make sure that multiple motifs are accepted as long as the total
        // score is increasing
        for (auto &result : results) {
          size_t task_idx = 0;  // TODO FIX choose the right task; for now
                                // assume it is the first one
          if (*result.training.state.scores[task_idx].rbegin()
              == best_score_for_this_motif) {
            HMM best_model = result.model;
            options = result.options;
            result.create_symlinks();
          }
        }
      }
    }
  }
  return hmm;
}

vector<HMM> cross_validation(const Data::Collection &all_data,
                             const Options::HMM &options) {
  vector<HMM> hmms;
  for (size_t cross_validation_iteration = 0;
       cross_validation_iteration < options.cross_validation_iterations;
       cross_validation_iteration++) {
    if (options.verbosity >= Verbosity::info
        and options.cross_validation_iterations > 1)
      cout << "Doing cross-validation " << (cross_validation_iteration + 1)
           << " of " << options.cross_validation_iterations << "." << endl;

    Options::HMM opt = options;

    if (options.cross_validation_freq < 1)
      opt.label += ".cv" + to_string(cross_validation_iteration);

    Data::Collection training_data, test_data;
    prepare_cross_validation(all_data, training_data, test_data,
                             options.cross_validation_freq, options.verbosity);
    HMM hmm = doit(all_data, training_data, test_data, opt);
    hmms.push_back(hmm);
  }
  return hmms;
}

void perform_analysis(Options::HMM &options) {
  if (options.verbosity >= Verbosity::verbose)
    cout << "Loading sequences." << endl;

  Data::Collection collection(options.paths, options.revcomp, options.n_seq);
  if (not options.dont_save_shuffle_sequences) {
    auto paths = collection.save_shuffle_sequences(options.label);
    if (options.verbosity >= Verbosity::info)
      for (auto &path : paths)
        cout << "Saving generated shuffle sequences to " << path << "." << endl;
  }

  check_data(collection, options);

  if (options.cross_validation_iterations == 0
      or options.cross_validation_freq == 1) {
    options.cross_validation_freq = 1;
    options.cross_validation_iterations = 1;
  }
  vector<HMM> hmms = cross_validation(collection, options);
}
