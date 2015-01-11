
#ifndef DREME_HPP

#include <cstddef>
#include <list>
#include <string>
#include <exception>

namespace Dreme {
std::list<std::pair<std::string, double>> parse_dreme_output(
    const std::string &dir);

std::list<std::pair<std::string, double>> run(
    const std::string &path1, const std::string &path2 = "",
    size_t min_size = 0, size_t max_size = 0, bool revcomp = false,
    size_t n_motifs = 0, bool remove_temp_dir = true);

namespace Exception {
struct BinaryNotFound : public std::exception {
  const char *what() const noexcept;
};

struct InvalidLengths : public std::exception {
  InvalidLengths(size_t s1, size_t s2);
  const char *what() const noexcept;
  size_t min_size, max_size;
};

struct ReturnValueNonZero : public std::exception {
  const char *what() const noexcept;
  size_t return_value;
};
}
}

#endif
