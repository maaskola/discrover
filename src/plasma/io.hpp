
#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>


namespace IO {
  namespace Exception {
    namespace File {
      struct Existence : public std::runtime_error { 
        const std::string msg = "File does not exist: ";
        Existence(const std::string& path_) : std::runtime_error(msg + path_), path(path_) { };
        std::string path;
      };
      struct Access: public std::runtime_error { 
        const std::string msg = "File access failed: ";
        Access(const std::string& path_) : std::runtime_error(msg + path_), path(path_) { };
        std::string path;
      };
    }
  }
}

template <typename X>
void parse_file(const std::string &path, X fnc) throw(IO::Exception::File::Existence, IO::Exception::File::Access) {
  if(not boost::filesystem::exists(path))
    throw(IO::Exception::File::Existence(path));

  bool use_gzip = path.substr(path.size() - 3, 3) == ".gz";
  bool use_bzip2 = path.substr(path.size() - 4, 4) == ".bz2";
  std::ios_base::openmode flags = std::ios_base::in;
  if(use_gzip or use_bzip2)
    flags |= std::ios_base::binary;

  std::ifstream file(path, flags);
  if(not file)
    throw(IO::Exception::File::Access(path));
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  if(use_gzip)
    in.push(boost::iostreams::gzip_decompressor());
  if(use_bzip2)
    in.push(boost::iostreams::bzip2_decompressor());
  in.push(file);

  fnc(in);
};

#endif

