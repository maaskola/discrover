#include <string>
#include <sstream>
#include "exception.hpp"

using namespace std;

namespace Exception {
namespace Parsing  {
namespace NumberList {
InvalidCharacter::InvalidCharacter(const string &spec_, size_t pos_) :
  runtime_error(""), spec(spec_), pos(pos_) { };
const char *InvalidCharacter::what() const noexcept {
  stringstream ss;
  ss << "Found invalid character '" << spec[pos] << "' at position " << pos << " of number list specification '" << spec << "'." << endl;
  ss << "Please note that the format for the list specification only allows digits, '-', and ','." << endl;
  return ss.str().c_str();
}

MultipleRanges::MultipleRanges(const string &group_) :
  runtime_error(""), group(group_) { };
const char *MultipleRanges::what() const noexcept {
  stringstream ss;
  ss << "List format error: only one '-' is allowed in any group." << endl
    << "The offending group is '" << group << "'." << endl;
  return ss.str().c_str();
}
};
};
};
