#include "random_seed.hpp"

size_t generate_rng_seed() {
  size_t seed = time(0) ^ getpid(); // XOR Unix time & process ID
  mix_with_os_entropy(seed);
  return(seed);
}
