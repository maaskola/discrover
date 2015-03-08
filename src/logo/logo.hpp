#ifndef LOGO_HPP
#define LOGO_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include "options.hpp"

namespace Logo {
using column_t = std::vector<double>;
using matrix_t = std::vector<column_t>;

/** Draw PDF or PNG motif logos.
 * The widths argument can be used to direct horizontal scaling of states.
 * It must either be empty, or have the same length as the number of columns in
 * the matrix.
 */
std::vector<std::string> draw_logo(const matrix_t &matrix,
                                   const std::string &path,
                                   const Logo::Options &options,
                                   const std::vector<double> &widths
                                   = std::vector<double>());
}

namespace Exception {
namespace Logo {
struct InvalidWidthArgument : public std::runtime_error {
  InvalidWidthArgument(size_t found, size_t expected);
};
}
}

#endif
