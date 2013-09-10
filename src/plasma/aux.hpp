


#ifndef AUX_HPP
#define AUX_HPP

#include <string>
#include <vector>

template <typename Iter> void range_tolower(Iter beg, Iter end) {
  for(Iter iter = beg; iter != end; ++iter) {
    *iter = std::tolower(*iter);
  }
}

std::string string_tolower(const std::string & str);

template <class T> T hibit_orig(T n) {
  n |= (n >>   1);
  n |= (n >>   2);
  n |= (n >>   4);
  n |= (n >>   8);
  n |= (n >>  16);
  n |= (n >>  32);
  return n - (n >> 1);
};

template <class T> T hibit(T n) {
  if(n == 0)
    return(0);
  size_t bits = 1;
  while((n = (n >> 4)) > 0)
    bits++;
  return bits;
};

/** Parse a comma separated list of ranges.
 * Example: "1,2,5-7,2,5,30"  will yield {1, 2, 5, 6, 7, 2, 5, 30}. */
std::vector<std::string> tokenize(const std::string &s, const std::string &delim);
std::vector<size_t> parse_list(const std::string &s);

std::string sha1hash(const std::string &s);

#endif   /* ----- #ifndef AUX_HPP  ----- */

