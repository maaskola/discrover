#include "random_seed.hpp"
#include <random>

using namespace std;
size_t generate_rng_seed() {
  uniform_int_distribution<size_t> udist;
  random_device rng;
  return udist(rng);
}
