
#ifndef DREME_HPP

#include <cstddef>
#include <list>
#include <string>
#include <exception>

namespace Dreme {
  class BinaryNotFoundException: public std::exception {
    virtual const char* what() const throw();
  };

  class InvalidLengthsException: public std::exception {
    private:
      size_t min_size, max_size;
    public:
      InvalidLengthsException(size_t s1, size_t s2);
      virtual const char* what() const throw();
  };

  std::list<std::pair<std::string,double>> parse_dreme_output(const std::string &dir);

  std::list<std::pair<std::string,double>> run(const std::string &path1,
      const std::string &path2="",
      size_t min_size=0,
      size_t max_size=0,
      bool revcomp=false,
      size_t n_motifs=0,
      bool remove_temp_dir=true);
};

#endif

