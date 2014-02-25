
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "../aux.hpp"
#include "hmm.hpp"
#include "report.hpp"
#include "../timer.hpp"
#include "../plasma/find.hpp"
#include "../stats_config.hpp"

using namespace std;

void check_data(const Data::Collection &collection, const hmm_options &options)
{
  if(options.verbosity >= Verbosity::info)
    for(auto &series: collection) {
      cout << "Data collection has " << series.set_size << " sequences with a total size of " << series.seq_size << " nucleotides." << endl;
      for(auto &set: series) {
        cout << set.path << " has " << set.set_size << " sequences with a total size of " << set.seq_size << " nucleotides." << endl;
        if(options.verbosity >= Verbosity::verbose) {
          cout << "The SHA1 of this set is " << set.sha1 << endl;
          cout << "The first 3 sequences are:" << endl;
          size_t i = 0;
          for(auto &seq: set) {
            cout << ">" << seq.definition << endl << seq.sequence << endl;
            if(++i >= 3) break;
          }
        }
        if(options.verbosity >= Verbosity::debug)
          for(auto &seq: set)
            cout << ">" << seq.definition << endl << seq.sequence << endl;
      }
    }

// TODO re-enable warning
//  if(options.training_method != Training::Method::none and not is_generative(options.objective.measure) and options.objective.measure != Measure::rank_information) {
//    for(auto &series: collection) {
//      if(series.sets.size() < 2) {
//        cout << "Please note that for discriminative training you need to specify multiple sequence sets using -f." << endl;
//        exit(-1);
//      }
//    }
//  }
}


void train_evaluate(HMM &hmm, const Data::Collection &all_data, const Data::Collection &training_data, const Data::Collection &test_data, const hmm_options &options)
{
  // Define the training tasks
  Training::Tasks tasks = hmm.define_training_tasks(options);

  if(not tasks.empty())
    hmm.train(training_data, tasks, options);
  if(test_data.set_size != 0) {
    evaluate_hmm(hmm, training_data, "training", tasks, options);
    evaluate_hmm(hmm, test_data, "Test", tasks, options);
  }
  evaluate_hmm(hmm, all_data, "", tasks, options);
}

double get_expected_seq_size(const Data::Collection &collection)
{
  return(1.0 * collection.seq_size / collection.set_size);
}

string padding(size_t n, char c='n') {
  string s;
  for(size_t i = 0; i < n; i++)
    s += c;
  return(s);
}

/* FIXME remove
string pad(const string &s, size_t l, size_t r, char c='n') {
  return(padding(l,c) + s + padding(r,c));
}  */

vector<string> generate_wiggle_variants(const string &s,  size_t n, Verbosity verbosity) {
  if(verbosity >= Verbosity::debug)
    cout << "generate_wiggle_variants(" << s << ", " << n << ")" << endl;
  size_t l = s.size();
  vector<string> variants;
  variants.push_back(s);
  for(size_t i = 1; i <= min<size_t>(l, n); i++) {
    variants.push_back(padding(i) + s.substr(0, l - i));
    variants.push_back(s.substr(i, l - i) + padding(i));
  }
  if(verbosity >= Verbosity::debug)
    for(auto &x: variants)
      cout << "Wiggle variant: " << x << endl;
  return(variants);
}


HMM doit(const Data::Collection &all_data, const Data::Collection &training_data, const Data::Collection &test_data, const hmm_options &options_)
{
  hmm_options options = options_;

  if(options.verbosity >= Verbosity::debug)
    cout << "About to construct HMM." << endl;

  // initialize HMM
  HMM hmm(options.bg_order, options.verbosity, options.contingency_pseudo_count);

  // load HMM from file if specified
  size_t n_loaded = 0;
  for(auto &load_path: options.load_paths)
    if(n_loaded++ == 0)
      hmm = HMM(load_path, options.verbosity, options.contingency_pseudo_count);
    else
      hmm.add_motifs(HMM(load_path, options.verbosity, options.contingency_pseudo_count));

  if(options.verbosity >= Verbosity::verbose)
    cout << "Model after initialization = " << hmm << endl;

  hmm.switch_intermediate(options.store_intermediate);

  // train background
  if(n_loaded == 0 and not options.objectives.empty()) {
    hmm.initialize_bg_with_bw(training_data, options);

    if(options.verbosity >= Verbosity::verbose)
      cout << "Model after background learning = " << hmm << endl;
  }

  double expected_seq_size = get_expected_seq_size(all_data);
  if(options.verbosity >= Verbosity::info)
    cout << "The average sequence size is " << expected_seq_size << "nt." << endl;

  for(auto &spec: options.motif_specifications)
    switch(spec.kind) {
      case Specification::Motif::Kind::File:
        // load emission matrix and add it.
        hmm.add_motif(read_emission(spec.specification), expected_seq_size, options.lambda, spec.name, spec.insertions, options.left_padding, options.right_padding);
        break;
      case Specification::Motif::Kind::Seed:
        // use IUPAC regular expression
        hmm.add_motif(spec.specification, options.alpha, expected_seq_size, options.lambda, spec.name, spec.insertions, options.left_padding, options.right_padding);
        break;
      case Specification::Motif::Kind::Plasma:
        options.seeding.motif_specifications.push_back(spec);
        break;
    }

  if(find_if(begin(options.motif_specifications), end(options.motif_specifications), [](const Specification::Motif &motif) {
        return(motif.kind == Specification::Motif::Kind::Plasma);
        }) == end(options.motif_specifications)) {
    if(options.verbosity >= Verbosity::info)
      cout << "No automatic seeds are used." << endl;
    train_evaluate(hmm, all_data, training_data, test_data, options);
  } else {
    if(options.verbosity >= Verbosity::info)
      cout << "Determining seeds automatically." << endl;

    Seeding::DataCollection collection(training_data);

    if(options.verbosity >= Verbosity::debug) {
      for(auto &ser: training_data)
        for(auto &set: ser) {
          cerr << "HMM::Series " << ser.name << " set -> motifs:";
          for(auto &m: set.motifs)
            cerr << " " << m;
          cerr << endl;
        }
    }
    if(options.verbosity >= Verbosity::debug) {
      for(auto &ser: collection)
        for(auto &set: ser) {
          cerr << "Seeding::Series " << ser.name << " set -> motifs:";
          for(auto &m: set.motifs)
            cerr << " " << m;
          cerr << endl;
        }
    }

    Seeding::Plasma plasma(collection, options.seeding);

    if(options.verbosity >= Verbosity::debug)
      cout << "motif_specs.size() = " << plasma.options.motif_specifications.size() << endl;
    size_t plasma_motif_idx = 0;

    bool first_motif = true;

    // while there are motif specifications left
    while(plasma_motif_idx < plasma.options.motif_specifications.size()) {

      if(hmm.get_nmotifs() > 0) {
        cout << "Masking plasma data collection." << endl;
        plasma.collection.mask(hmm.compute_mask(training_data));
      }

      // consider the next motif specification
      auto motif_spec = options.motif_specifications[plasma_motif_idx];

      // determine the matching Plasma objective
      auto objectives = Training::corresponding_objectives(options.objectives, options.use_mi_to_seed);
      plasma.options.objectives = objectives;
      if(options.verbosity >= Verbosity::info) {
        cout << "Determining seed " << motif_spec << " with objective functions";
        for(auto &obj: objectives)
          cout << " " << obj;
        cout << "." << endl;
      }

      // find the next set of seeds with Plasma
      Seeding::Results plasma_results = plasma.find(motif_spec, objectives);

      // if no seeds are found, don't try to find more seeds
      if(plasma_results.empty())
        break;
      if(options.verbosity >= Verbosity::debug)
        cout << motif_spec.name << " result.size() = " << plasma_results.size() << endl;

      vector<pair<string, HMM>> learned_models;

      // seed and learn HMM parameters independently for each Plasma motif
      for(size_t seed_idx = 0; seed_idx < plasma_results.size(); seed_idx++) {
        string motif = plasma_results[seed_idx].motif;
        string name = motif_spec.name;

        cout << "Considering candidate motif " << name << ":" << motif << " and training to determine the HMM score." << endl;

        plasma_results[seed_idx].score = -numeric_limits<double>::infinity();

        // consider all wiggle variants
        for(auto variant: generate_wiggle_variants(motif, options.wiggle, options.verbosity)) {
          HMM model(hmm);
          if(options.verbosity >= Verbosity::info and options.wiggle > 0)
            cout << "Considering wiggle variant " << variant << " of candidate motif " << name << ":" << motif << " and training to determine the HMM score." << endl;
          model.add_motif(variant, options.alpha, expected_seq_size, options.lambda, name, plasma.options.motif_specifications[plasma_motif_idx].insertions, options.left_padding, options.right_padding);

          hmm_options options_(options);
          if(options_.long_names)
            options_.label += "." + variant;
          train_evaluate(model, all_data, training_data, test_data, options_);

          learned_models.push_back(make_pair(variant, model));
        }
      }

      if(not options.accept_multiple) {
        if(learned_models.size() == 1)
          hmm = learned_models.begin()->second;
        else {
          if(options.verbosity >= Verbosity::info)
            cout << "Evaluating learned models.";
          auto best_model = hmm;
          string best_seed = "";
          double best_score = -numeric_limits<double>::infinity();
          for(auto &learned: learned_models) {
            string seed = learned.first;
            auto model = learned.second;
            Training::Tasks tasks = model.define_training_tasks(options);
            // TODO consider using compute_score instead of compute_score_all_motifs
            double score = model.compute_score_all_motifs(training_data, tasks.begin()->measure, options.weighting);
            if(score > best_score) {
              best_score = score;
              best_seed = seed;
              best_model = model;
            }
          }
          if(options.verbosity >= Verbosity::info)
            cout << "Accepting seed " << best_seed << " with score " << best_score << endl;
          hmm = best_model;
        }
      } else {
        auto best_score = -numeric_limits<double>::infinity();
        bool ok = true;

        cout << "Non-augmented model: " << hmm << endl;
        // Training::Tasks tasks = hmm.define_training_tasks(options);
        // double previous_score = hmm.compute_score(masked_training_data, *tasks.begin(), options.weighting);
        while(ok and not learned_models.empty()) {
          auto masked_training_data = training_data;

          // mode == 0:
          //   do not mask occurrences of previously identified motifs;
          //   add candidate motif to current model
          //   score composite model consisting of candidate motif and previously identified motifs on unmasked data
          // mode == 1:
          //   mask occurrences of previously identified motifs;
          //   score candidate motif as single motif model on masked data
          // mode == 2:
          //   mask occurrences of previously identified motifs;
          //   add candidate motif to current model
          //   score composite model consisting of candidate motif and previously identified motifs on masked data
          const size_t mode = 2;
          const bool relearn_before_eval = true;

          if(mode == 1 or mode == 2) {
            best_score = -numeric_limits<double>::infinity();
            if(not first_motif) {
              cout << "Masking earlier motifs" << endl;
              masked_training_data.mask(hmm.compute_mask(masked_training_data));
            }
          }

          size_t best_index = 0;
          auto best_model = hmm;
          string best_seed = "";
          bool updated = false;

          size_t index = 0;

          for(auto &learned: learned_models) {
            string seed = learned.first;
            auto learned_model = learned.second;

            auto model = learned_model;
            if(mode == 0 or mode == 2) {
              model = hmm;

              // TODO add the motif from the learned model to model
              model.add_motifs(learned_model, false);
              // model.add_motifs(learned_model, true);

              // TODO adapt transition probabilities
              // TODO or do complete relearning
              // if(relearn_before_eval)
            } else
              model = learned_model;


            if(false) { // currently not needed; might move this code done to the relearning part
              // TODO make sure the task is actually referring to the newly added motif so that is the motif that is scored
              Training::Tasks tasks = model.define_training_tasks(options);
              // TODO this means we want to re-write the discriminative objectives' targets so that they refer only to the last = most recently added motif
              auto group_it = hmm.groups.rbegin();
              auto not_in_group = [&group_it](size_t v) { return(find(begin(group_it->states), end(group_it->states), v) == end(group_it->states)); };
              for(auto &task: tasks)
                if(Measures::is_discriminative(task.measure)) {
                  remove_if(task.targets.emission.begin(), task.targets.emission.end(), not_in_group);
                  remove_if(task.targets.transition.begin(), task.targets.transition.end(), not_in_group);
                }
            }

            // potentially: regardless of the objective function chosen for training, one might use the MICO p-value for selection!
            const bool use_mico_pvalue = true;
            double score;
            std::vector<size_t> groups_to_score;
            groups_to_score.push_back(model.groups.size() - 1); // only add the most recently added group
            if(use_mico_pvalue)
              score = model.compute_score(masked_training_data, Measures::Continuous::Measure::MutualInformation, options.weighting, groups_to_score);
            else
              // TODO actually select the specific task - above we do fix all conceivable tasks...
              score = model.compute_score(masked_training_data, model.define_training_tasks(options).begin()->measure, options.weighting, groups_to_score);
            if(use_mico_pvalue) {
              cout << "mi = " << score << endl;
              double n = masked_training_data.set_size; // TODO fix this with regards to pseudo counts and exact reference to the relevant contrast
              double df = 1;
              size_t motif_len = seed.size();
              double g = calc_g_test_from_mi(score, n);
              cout << "g = " << g << endl;
              double log_p = pchisq(g, df, false, true);
              cout << "log p(g) = " << log_p << endl;
              double cor_log_p = log(149) * motif_len + log_p;
              cout << "corrected log p(g) = " << cor_log_p << endl;
              score = - cor_log_p;
              cout << "score = " << score << endl;
            }

            cout << "Augmented model: " << model << endl;
            cout << "Score of the model augmented by " << seed << " has a score of " << score << endl;
            if(score > best_score) {
              if(not use_mico_pvalue or score > 0) {
                ok = true;
                updated = true;
                best_model = model;
                best_seed = seed;
                best_score = score;
                best_index = index;
              }
            }
            index++;
          }

          if(not updated) {
            cout << "We did not find an improved model." << endl;
            ok = false;
          } else {
            if(mode == 0 or mode == 2)
              hmm = best_model;

            if(options.verbosity >= Verbosity::info)
              cout << "Accepting seed " << best_seed << " with score " << best_score << endl;

            if(mode == 1)
              hmm.add_motifs(best_model, false);

            hmm_options options_(options);
            if(options.long_names)
              options.label += "." + best_seed;

            // TODO: remember that learning and evaluation is to a large degree based on groups - not motif names; thus there is some inefficiencies
            // TODO: learning should perhaps only adapt transitions? this would be very fast...
            // TODO:   otherwise: if learning should also adapt emissions, then perhaps only those of the new motif?
            // TODO:   in this case: it might be done on the masked sequences
            train_evaluate(hmm, all_data, training_data, test_data, options);

            learned_models.erase(begin(learned_models) + best_index);
            first_motif = false;
            ok = true;
          }
        }
      }

      plasma_motif_idx++;
    }
  }
  return(hmm);
}

vector<HMM> cross_validation(const Data::Collection &all_data, const hmm_options &options)
{
  vector<HMM> hmms;
  for(size_t cross_validation_iteration = 0; cross_validation_iteration < options.cross_validation_iterations; cross_validation_iteration++) {
    if(options.verbosity >= Verbosity::info and options.cross_validation_iterations > 1)
      cout << "Doing cross-validation " << (cross_validation_iteration + 1) << " of " << options.cross_validation_iterations << "." << endl;

    hmm_options opt = options;

    if(options.cross_validation_freq < 1)
      opt.label += ".cv" + boost::lexical_cast<string>(cross_validation_iteration);

    Data::Collection training_data, test_data;
    prepare_cross_validation(all_data, training_data, test_data, options.cross_validation_freq, options.verbosity);
    HMM hmm = doit(all_data, training_data, test_data, opt);
    hmms.push_back(hmm);
  }
  return(hmms);
}

void perform_analysis(hmm_options &options)
{
  if(options.verbosity >= Verbosity::verbose)
    cout << "Loading sequences." << endl;

  Data::Collection collection(options.paths, options.revcomp, options.n_seq);

  check_data(collection, options);

  if(options.cross_validation_iterations == 0 or options.cross_validation_freq == 1) {
    options.cross_validation_freq = 1;
    options.cross_validation_iterations = 1;
  }
  vector<HMM> hmms = cross_validation(collection, options);
}

