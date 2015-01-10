#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <stdexcept>

namespace Exception {
namespace Parsing  {
namespace NumberList {
struct InvalidCharacter : public std::runtime_error {
  InvalidCharacter(const std::string &spec, size_t pos);
  std::string spec;
  size_t pos;
  const char *what() const noexcept;
};
struct MultipleRanges : public std::runtime_error {
  MultipleRanges(const std::string &group);
  std::string group;
  const char *what() const noexcept;
};
};
};
namespace BitMask {
struct TooManyMotifs : public std::length_error {
  TooManyMotifs (size_t n_motifs, size_t max_motifs) const noexcept;
  size_t n_motifs;
  size_t max_motifs;
};
};
};

#endif
