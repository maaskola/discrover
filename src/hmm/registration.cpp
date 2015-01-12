
#include <iostream>
#include "registration.hpp"
#include "hmm.hpp"

using namespace std;

Registration::Registration(Verbosity v) : verbosity(v){};

void Registration::add_dataset(const Data::Set &dataset, double class_prior) {
  if (datasets.find(dataset.sha1) == end(datasets)) {
    if (verbosity >= Verbosity::verbose)
      cout << "register_data_set(dataset.path=" << dataset.path
           << (dataset.is_shuffle ? " shuffle" : " ")
           << ", sha1=" << dataset.sha1 << ", class_prior=" << class_prior
           << ")" << endl;
    Sample sample = {dataset, class_prior, unordered_map<bitmask_t, double>()};
    datasets[dataset.sha1] = sample;
  } else if (verbosity >= Verbosity::verbose)
    cout << "Not registering already present data_set(dataset.path="
         << dataset.path << (dataset.is_shuffle ? " shuffle" : " ")
         << ", sha1=" << dataset.sha1 << ", class_prior=" << class_prior << ")"
         << endl;
}

void Registration::add_bitmask(const std::string name, bitmask_t present,
                               double motif_p1, double motif_p2) {
  for (auto &sample : datasets) {
    if (sample.second.motif_prior.find(present)
        == end(sample.second.motif_prior)) {
      if (sample.second.spec.motifs.find(name)
          != end(sample.second.spec.motifs))
        sample.second.motif_prior[present] = motif_p1;
      else
        sample.second.motif_prior[present] = motif_p2;
      if (verbosity >= Verbosity::verbose)
        cout << "Adding motif prior for motif " << name << " / " << present
             << " in data set " << sample.second.spec.path << " = "
             << sample.second.motif_prior[present] << endl;
    } else if (verbosity >= Verbosity::verbose)
      cout << "Not adding motif prior for motif " << name << " / " << present
           << " in data set " << sample.second.spec.path << " = "
           << sample.second.motif_prior[present] << endl;
  }
}

double Registration::Sample::get_motif_prior(bitmask_t present) const {
  auto x = motif_prior.find(present);
  if (x == end(motif_prior))
    throw Exception::Registration::UnregisteredMotifGroup(present);
  else
    return x->second;
}

double Registration::compute_marginal_motif_prior(bitmask_t present) const {
  if (verbosity >= Verbosity::debug)
    cout << "compute_marginal_motif_prior(present=" << present << ")" << endl;
  double marginal_motif_prior = 0;
  for (auto &x : datasets)
    marginal_motif_prior += x.second.class_prior
                            * x.second.get_motif_prior(present);
  return marginal_motif_prior;
}

double Registration::get_class_motif_prior(const string &sha1,
                                           bitmask_t present) const {
  if (verbosity >= Verbosity::debug)
    cout << "get_class_motif_prior(sha1=" << sha1 << ", present=" << present
         << ")" << endl;

  auto cparms = datasets.find(sha1);
  if (cparms == end(datasets))
    throw Exception::Registration::UnregisteredDataSet(sha1);
  auto x = cparms->second.motif_prior.find(present);
  if (x == end(cparms->second.motif_prior))
    throw Exception::Registration::UnregisteredMotifGroup(present);
  if (verbosity >= Verbosity::debug)
    cout << "get_class_motif_prior(sha1=" << sha1 << ", present=" << present
         << ") = " << x->second << endl;
  return x->second;
}

double Registration::get_class_prior(const string &sha1) const {
  if (verbosity >= Verbosity::debug)
    cout << "get_class_prior(sha1=" << sha1 << ")" << endl;

  auto cparms = datasets.find(sha1);
  if (cparms == end(datasets))
    throw Exception::Registration::UnregisteredDataSet(sha1);
  double x = cparms->second.class_prior;
  if (verbosity >= Verbosity::debug)
    cout << "get_class_prior(sha1=" << sha1 << ") = " << x << endl;
  return x;
}

namespace Exception {
namespace Registration {
UnregisteredMotifGroup::UnregisteredMotifGroup(const bitmask_t &present_)
    : exception(), present(present_) {};
const char *UnregisteredMotifGroup::what() const noexcept {
  string msg = "Error: could not find class parameters for motif group:"
               + present.to_string() + ".";
  return msg.c_str();
}
UnregisteredDataSet::UnregisteredDataSet(const string &sha1_)
    : exception(), sha1(sha1_) {};
const char *UnregisteredDataSet::what() const noexcept {
  string msg = "Error: could not find class parameters of sequences with sha1 "
               + sha1 + ".";
  return msg.c_str();
}
}
}
