#include <string>
#include <vector>
#include "options.hpp"

namespace Logo {
using column_t = std::vector<double>;
using matrix_t = std::vector<column_t>;

std::vector<std::string> draw_logo(const matrix_t &matrix, const std::string &path, const Logo::Options &options);
};
