
#include "randomization_test.hpp"
#include "score.hpp"

using namespace std;
typedef Seeding::Stats::OccurrenceCounts counts_t;

counts_t generate_sample(const Seeding::DataCollection &collection, const counts_t &counts, Verbosity verbosity) {
  counts_t sample(counts.size(), 0);

  if(verbosity >= Verbosity::debug)
    cout << "Sample generation"
      << "original counts = " << counts << endl
      << "sampled counts = " << sample << endl;

  size_t start_idx = 0;

  for(auto &series: collection) {
    size_t stop_idx = start_idx + series.sets.size();

    if(verbosity >= Verbosity::debug)
      cout << "Considering series with count start index " << start_idx << " and stop index " << stop_idx << endl;

    size_t successes = 0;
    for(size_t i = start_idx; i < stop_idx; i++)
      successes += counts(i);

    size_t remaining = series.set_size;

    if(verbosity >= Verbosity::debug)
      cout << "In total there are " << successes << " in this series of " << remaining << " sequences." << endl;

    size_t k = 0;
    for(auto &set: series) {
      for(auto &seq: set) {
        double r = 1.0 * rand() / RAND_MAX;
        if(r < 1.0 * successes / remaining) {
          successes--;
          sample[k]++;
        }
        remaining--;
      }
      k++;
    }

    assert(successes == 0);
    assert(remaining == 0);

    start_idx = stop_idx;
  }

  return(sample);
}

bool randomization_test(const Seeding::DataCollection &collection, const counts_t &counts, size_t nr_tests, double score, const Seeding::Options &options, const Seeding::Objective &objective, size_t length, size_t degeneracy) {
  if(options.verbosity >= Verbosity::verbose)
    cout << "Randomization test" << endl;
  for(size_t i = 0; i < nr_tests; i++) {

    counts_t sample = generate_sample(collection, counts, options.verbosity);

    double sample_score = compute_score(collection, sample, options, objective, length, degeneracy);

    bool success = sample_score < score;

    if(options.verbosity >= Verbosity::debug)
      cout << "original counts = " << counts << endl
        << "sampled counts = " << sample << endl
        << "original score = " << score << endl
        << "sampled score = " << sample_score << endl
        << "Success = " << success << endl;

    if(not success) {
      if(options.verbosity >= Verbosity::verbose)
        cout << "Failed after " << (i+1) << " tests." << endl;

      return(false);
    }
  }

  if(options.verbosity >= Verbosity::verbose)
    cout << "Randomization test succeeded " << nr_tests << " times." << endl;

  return(true);
}
