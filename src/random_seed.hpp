#ifndef RANDOM_SEED_HPP

#include <ctime>          // for time()
#include <fstream>
#include <boost/filesystem.hpp>

template <typename X> X entropy_from_os(const std::string &source="/dev/urandom") {
  std::ifstream f(source.c_str());
  X myRandom;
  f.read(reinterpret_cast<char*>(&myRandom), sizeof(myRandom));
  return(myRandom);
}

template <typename X> void mix_with_os_entropy(X &x, const std::string &source="/dev/urandom") {
  if(boost::filesystem::exists(source)) {
    X y = entropy_from_os<X>(source);
    x = x ^ y;
  }
}

size_t generate_rng_seed();

#endif
