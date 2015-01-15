#include "random_seed.hpp"
#include "random_distributions.hpp"

using namespace std;
size_t generate_rng_seed() {
  random_device rng;
  return RandomDistribution::Uniform(rng);
}
