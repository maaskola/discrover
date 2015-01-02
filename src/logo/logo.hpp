#include <string>
#include <vector>

namespace Logo {
using column_t = std::vector<double>;
using matrix_t = std::vector<column_t>;

enum class output_t {
  PDF,
  PNG
};

bool draw_logo(const matrix_t &matrix, const std::string &path, output_t kind);
};
