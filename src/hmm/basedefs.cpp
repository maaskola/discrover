
#include <fstream>
#include <set>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "../plasma/fasta.hpp"
#include "aux.hpp"
#include "basedefs.hpp"

using namespace std;

vector<string> get_paths(const std::vector<Specification::DataSet> &specs)
{
  vector<string> paths;
  for(auto &s: specs)
    paths.push_back(s.path);
  return(paths);
}

void prepare_cross_validation(const Data::Series &data_sets, Data::Series &training_data, Data::Series &test_data, double cross_validation_freq, Verbosity verbosity)
{
  if(cross_validation_freq == 1)
    training_data = data_sets;
  else {
    for(auto &data: data_sets) {
      if(verbosity >= Verbosity::verbose)
        std::cerr << "Splitting " << data.path << " into training and test data." << std::endl;
      Data::Set training, test;

      training.sha1 = data.sha1;
      test.sha1 = data.sha1;

      training.path = data.path;
      test.path = data.path;

      training.motifs = data.motifs;
      test.motifs = data.motifs;

      training.series = data.series;
      test.series = data.series;

      for(auto &seq: data) {
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
        std::cerr << "Training data set size of " << data.path << " = " << training.set_size << std::endl;
        std::cerr << "Test data set size of " << data.path << " = " << test.set_size << std::endl;
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

void prepare_cross_validation(const Data::Collection &data_sets, Data::Collection &training_data, Data::Collection &test_data, double cross_validation_freq, Verbosity verbosity)
{
  for(auto &series: data_sets) {
    Data::Series training, test;
    training.name = series.name;
    test.name = series.name;

    prepare_cross_validation(series, training, test, cross_validation_freq, verbosity);

    training_data.series.push_back(training);
    training_data.seq_size += training.seq_size;
    training_data.set_size += training.set_size;

    test_data.series.push_back(test);
    test_data.seq_size += test.seq_size;
    test_data.set_size += test.set_size;
  }
}


namespace Training {
  State::State(size_t n) : center(-9), scores(n)
  {
  }
}

