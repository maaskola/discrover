#include "random_seed.hpp"
#include <random>

using namespace std;
size_t generate_rng_seed() {
  uniform_int_distribution<size_t> r_unif;
  random_device rng;
  return r_unif(rng);
}
