#include <string>

namespace logo {

enum class output_t {
  PDF,
  PNG
};

bool draw_logo(const std::string &path, output_t kind);
};
