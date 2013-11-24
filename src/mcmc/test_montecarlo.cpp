#include <iostream>
#include <cstddef>
#include "montecarlo.hpp"

using namespace std;

int main(int argc, const char **argv) {
//  size_t n_samples = 1000000;
  size_t n_samples = 1000;
  if(argc == 2)
    n_samples = atoi(argv[1]);

//  MCMC::Generator<double> gen;
//  MCMC::Evaluator<double> eval;
  MCMC::MonteCarlo<double> mcmc(Verbosity::info);
  /*
  list<double> res = mcmc.run(1.0, 50, n_samples);
  size_t idx = 0;
  for(auto x : res)
    cout << idx++ << "\t" << x << endl;
    */
  vector<double> t = {10, 5, 1, 0.5, 0.1};
  vector<double> init = {50, 50, 50, 50, 50};
  auto res2 = mcmc.parallel_tempering(t, init, n_samples);
  size_t k = 0;
  for(auto x : res2) {
    size_t idx = 0;
    for(auto y : x)
      cout << k << "\t" << idx++ << "\t" << y.first << "\t" << y.second << endl;
    k++;
  }
  return(EXIT_SUCCESS);
}

