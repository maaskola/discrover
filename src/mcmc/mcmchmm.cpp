
#ifndef MCMCHMM_HPP
#define MCMCHMM_HPP

#include "../hmm/hmm.hpp"
#include "montecarlo.hpp"

namespace MCMC {
template <size_t alphabet_size>
class Generator<HMM<alphabet_size>> {
private:
  TrainingTargets train_what;
  size_t n_indels;
  size_t min_size;
  size_t max_size;

public:
  Generator(TrainingTargets what, size_t indels, size_t min_s, size_t max_s)
      : train_what(what),
        n_indels(indels),
        min_size(min_size),
        max_size(max_size){};
  HMM<alphabet_size> generate(const HMM<alphabet_size> &hmm) const {
    return hmm.random_variant(what, indels, min_size, max_size);
  };
};
template <size_t alphabet_size>
class Evaluator<HMM<alphabet_size>> {
private:
  std::vector<DataSet> &signal_data;
  std::vector<DataSet> &control_data;
  ObjectiveFunction objFunc;
  size_t signal_data_size;
  size_t control_data_size;
  std::vector<size_t> signal_data_sizes;
  std::vector<size_t> control_data_sizes;

public:
  Evaluator(const std::vector<DataSet> &signal,
            const std::vector<DataSet> &control, ObjectiveFunction fnc)
      : signal_data(signal),
        control_data(control),
        objFund(fnc),
        signal_data_size(0),
        control_data_size(0),
        signal_data_sizes(),
        control_data_sizes() {
    signal_data_sizes = get_dataset_sizes(signal, fnc);
    control_data_sizes = get_dataset_sizes(control, fnc);
    for (auto s : signal_data_sizes)
      signal_data_size += s;
    for (auto s : control_data_sizes)
      control_data_size += s;
  };

  double evaluate(const HMM<alphabet_size> &hmm) const {
    double score = -hmm.compute_score(signal_data, control_data, objFunc,
                                      signal_data_size, control_data_size);
    if (verbosity >= Verbosity::info)
      std::cout << "initial score = " << score << std::endl;
    return score;
  };
}
}

#endif
