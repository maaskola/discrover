
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

double train_hmm(HMM &hmm,
    const Data::Collection &training_data,
    const Training::Tasks &tasks,
    const hmm_options &options)
{
  if(options.verbosity >= Verbosity::verbose)
    cout << "Model to be evaluated = " << hmm << endl;
  double delta = 0;
  if(tasks.empty()) {
    if(options.verbosity >= Verbosity::info)
      cout << "Not performing training because no training tasks were specified." << endl;
  } else {
    bool any_found = false;
    for(auto &group: hmm.groups)
      for(auto &task: tasks)
        if(task.motif_name == group.name) {
          any_found = true;
        }
    if(not any_found) {
      if(options.verbosity >= Verbosity::info)
        cout << "Skipping training because no motifs specified in the tasks have corresponding states in the HMM." << endl;
    } else {
      if(options.verbosity >= Verbosity::info)
        cout << "Performing training." << endl;
      Timer learning_timer;

      if(options.verbosity >= Verbosity::verbose)
        cout << "Registering data sets for class based HMMs." << endl;

      for(auto &series: training_data)
        for(auto &data_set: series)
          hmm.register_dataset(data_set, (1.0*data_set.set_size)/training_data.set_size, options.conditional_motif_prior1, options.conditional_motif_prior2);

      delta = hmm.train(training_data, tasks, options);
      if(options.verbosity >= Verbosity::verbose)
        cout << endl << "The parameters changed by an L1-norm of " << delta << endl;

      double time = learning_timer.tock();
      if(options.timing_information)
        cerr << "Learning: " << time << " micro-seconds" << endl;

      if(options.verbosity >= Verbosity::debug)
        cout << "HMM after training:" << endl
          << hmm << endl;

      string store_to = options.label + ".hmm";

      if(options.verbosity >= Verbosity::info)
        cout << endl << "Parameters stored in " << store_to << endl;

      ofstream os(store_to.c_str());
      hmm.serialize(os, options.exec_info);
    }
  }
  return(delta);
}

void initialize_bg_with_bw(HMM &hmm, const Data::Collection &collection, const hmm_options &options)
{
  if(options.verbosity >= Verbosity::info)
    cout << "Initializing background of order " << options.bg_order << " with Baum-Welch algorithm." << endl;

  hmm_options bg_options = options;
  if(options.verbosity == Verbosity::info)
    bg_options.verbosity = Verbosity::error;

  Timer timer;
  hmm.train_background(collection, bg_options);
  double time = timer.tock();

  if(options.timing_information)
    cerr << "Background learning: " << time << " micro-seconds" << endl;
}

 
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
    train_hmm(hmm, training_data, tasks, options);
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
    std::cout << "generate_wiggle_variants(" << s << ", " << n << ")" << std::endl;
  size_t l = s.size();
  vector<string> variants;
  variants.push_back(s);
  for(size_t i = 1; i <= min<size_t>(l, n); i++) {
    variants.push_back(padding(i) + s.substr(0, l - i));
    variants.push_back(s.substr(i, l - i) + padding(i));
  }
  if(verbosity >= Verbosity::debug)
    for(auto &x: variants)
      std::cout << "Wiggle variant: " << x << std::endl;
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
    initialize_bg_with_bw(hmm, training_data, options);

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

  if(std::find_if(begin(options.motif_specifications), end(options.motif_specifications), [](const Specification::Motif &motif) {
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
    auto lesser_score = [](const Seeding::Result &a, const Seeding::Result &b) { return(a.score < b.score); };

    while(not plasma.options.motif_specifications.empty()) {
      if(options.verbosity >= Verbosity::debug)
        cout << "motif_specs.size() = " << plasma.options.motif_specifications.size() << endl;
      if(hmm.get_nmotifs() > 0)
        plasma.collection.mask(hmm.compute_mask(training_data));
      Seeding::Results all_plasma_results;
      size_t plasma_motif_idx = 0;
      while(plasma_motif_idx < plasma.options.motif_specifications.size()) {
        auto motif_spec = options.motif_specifications[plasma_motif_idx];
        auto objectives = Training::corresponding_objectives(options.objectives, options.use_mi_to_seed);
        plasma.options.objectives = objectives;
        if(options.verbosity >= Verbosity::info) {
          cout << "Determining seed " << motif_spec << " with objective functions";
          for(auto &obj: objectives)
            cout << " " << obj;
          cout << "." << endl;
        }
        Seeding::Results plasma_results = plasma.find(motif_spec, objectives);
        if(plasma_results.empty())
          break;
        if(options.verbosity >= Verbosity::debug)
          cout << motif_spec.name << " result.size() = " << plasma_results.size() << endl;

        if(options.model_choice == ModelChoice::HMMScore) {
          for(size_t seed_idx = 0; seed_idx < plasma_results.size(); seed_idx++) {
            string motif = plasma_results[seed_idx].motif;
            string name = motif_spec.name;

            std::cout << "Considering candidate motif " << name << ":" << motif << " and training to determine the HMM score." << std::endl;

            plasma_results[seed_idx].score = -std::numeric_limits<double>::infinity();

            for(auto variant: generate_wiggle_variants(motif, options.wiggle, options.verbosity)) {
              HMM hmm_(hmm);
              hmm_options options_(options);
              if(options_.verbosity >= Verbosity::info and options_.wiggle > 0)
                std::cout << "Considering wiggle variant " << variant << " of candidate motif " << name << ":" << motif << " and training to determine the HMM score." << std::endl;
              hmm_.add_motif(variant, options_.alpha, expected_seq_size, options_.lambda, name, plasma.options.motif_specifications[plasma_motif_idx].insertions, options.left_padding, options.right_padding);

              if(options_.long_names)
                options_.label += "." + variant;
              Training::Tasks tasks = hmm_.define_training_tasks(options_);
              train_evaluate(hmm_, all_data, training_data, test_data, options_);
              double score = hmm_.compute_score(training_data, *tasks.begin(), options_.weighting);
              if(score > plasma_results[seed_idx].score) {
                plasma_results[seed_idx].motif = variant;
                plasma_results[seed_idx].score = score;
              }
            }
          }
        }

        all_plasma_results.push_back(*max_element(plasma_results.begin(), plasma_results.end(), lesser_score));
        if(options.verbosity >= Verbosity::debug)
          cout << "all.size() = " << all_plasma_results.size() << endl;
        plasma_motif_idx++;
      }
      if(options.verbosity >= Verbosity::debug)
        cout << "found all" << endl;
      if(all_plasma_results.empty()) {
        if(options.verbosity >= Verbosity::info)
          cout << "Unable to find any seeds." << endl;
        plasma.options.motif_specifications.clear();
      } else {
        auto best_iter = max_element(all_plasma_results.begin(), all_plasma_results.end(), lesser_score);
        size_t best_idx = best_iter - all_plasma_results.begin();

        if(options.verbosity >= Verbosity::debug)
          cout << "best_idx = " << best_idx << endl;
        string motif = best_iter->motif;
        string name = plasma.options.motif_specifications[best_idx].name;

        if(options.verbosity >= Verbosity::info)
          cout << "Accepting seed " << name << ":" << motif << endl;

        hmm.add_motif(motif, options.alpha, expected_seq_size, options.lambda, name, plasma.options.motif_specifications[best_idx].insertions, options.left_padding, options.right_padding);
        if(options.long_names)
          options.label += "." + motif;

        plasma.options.motif_specifications.erase(plasma.options.motif_specifications.begin() + best_idx);
      }
      if(options.simultaneity == Training::Simultaneity::Sequential)
        train_evaluate(hmm, all_data, training_data, test_data, options);
    }
    if(options.simultaneity == Training::Simultaneity::Simultaneous and (options.model_choice != ModelChoice::HMMScore or not plasma.options.motif_specifications.empty()))
      train_evaluate(hmm, all_data, training_data, test_data, options);
  }
  return(hmm);
}

std::vector<HMM> cross_validation(const Data::Collection &all_data, const hmm_options &options)
{
  std::vector<HMM> hmms;
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
  std::vector<HMM> hmms = cross_validation(collection, options);
}

