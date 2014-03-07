
#include <fstream>
#include <set>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "../plasma/fasta.hpp"
#include "../aux.hpp"
#include "basedefs.hpp"

using namespace std;

vector<string> get_paths(const vector<Specification::Set> &specs)
{
  vector<string> paths;
  for(auto &s: specs)
    paths.push_back(s.path);
  return(paths);
}

void prepare_cross_validation(const Data::Contrast &contrast, Data::Contrast &training_data, Data::Contrast &test_data, double cross_validation_freq, Verbosity verbosity)
{
  if(cross_validation_freq == 1)
    training_data = contrast;
  else {
    for(auto &dataset: contrast) {
      if(verbosity >= Verbosity::verbose)
        cerr << "Splitting " << dataset.path << " into training and test data." << endl;
      Data::Set training, test;

      training.sha1 = dataset.sha1;
      test.sha1 = dataset.sha1;

      training.path = dataset.path;
      test.path = dataset.path;

      training.motifs = dataset.motifs;
      test.motifs = dataset.motifs;

      training.contrast = dataset.contrast;
      test.contrast = dataset.contrast;

      for(auto &seq: dataset) {
        double p = 1.0 * rand() / RAND_MAX;
        if(p >= cross_validation_freq) {
          test.sequences.push_back(seq);
          test.set_size += 1;
          test.seq_size += seq.sequence.size();
        } else {
          training.sequences.push_back(seq);
          training.set_size += 1;
          training.seq_size += seq.sequence.size();
        }
      }
      if(verbosity >= Verbosity::verbose) {
        cerr << "Training data set size of " << dataset.path << " = " << training.set_size << endl;
        cerr << "Test data set size of " << dataset.path << " = " << test.set_size << endl;
      }

      training_data.sets.push_back(training);
      training_data.seq_size += training.seq_size;
      training_data.set_size += training.set_size;

      test_data.sets.push_back(test);
      test_data.seq_size += test.seq_size;
      test_data.set_size += test.set_size;
    }
  }
}

void prepare_cross_validation(const Data::Collection &collection, Data::Collection &training_data, Data::Collection &test_data, double cross_validation_freq, Verbosity verbosity)
{
  for(auto &contrast: collection) {
    Data::Contrast training, test;
    training.name = contrast.name;
    test.name = contrast.name;

    prepare_cross_validation(contrast, training, test, cross_validation_freq, verbosity);

    training_data.contrasts.push_back(training);
    training_data.seq_size += training.seq_size;
    training_data.set_size += training.set_size;

    test_data.contrasts.push_back(test);
    test_data.seq_size += test.seq_size;
    test_data.set_size += test.set_size;
  }
}


namespace Training {
  State::State(size_t n) : center(-9), scores(n)
  {
  }
}

