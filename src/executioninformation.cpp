
#include <cstdio>
#include <ctime>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "executioninformation.hpp"

using namespace std;

string cmdline(int argc, const char **argv) {
  string cmd;
  for (int i = 0; i < argc; i++)
    cmd += (i != 0 ? " " : "") + string(argv[i]);
  return cmd;
}

ExecutionInformation generate_exec_info(const string &name,
                                        const string &hmm_version,
                                        const string &cmdline) {
  time_t rawtime;
  time(&rawtime);

  string datetime = ctime(&rawtime);
  datetime = datetime.substr(0, datetime.size() - 1);

  string dir = boost::filesystem::initial_path().string();

  string program_name = boost::filesystem::path(name).filename().string();

  ExecutionInformation exec_info
      = {program_name, hmm_version, cmdline, datetime, dir};
  return exec_info;
}

string generate_random_label(const string &prefix, size_t n_rnd_char,
                             Verbosity verbosity) {
  random_device rng;
  uniform_int_distribution<char> r_char('a', 'z');
  using namespace boost::posix_time;

  string label = prefix;

  // TODO perhaps add the process ID -> getpid()
  if (verbosity >= Verbosity::debug)
    cout << "Generating random label with prefix " << prefix << " and "
         << n_rnd_char << " random characters." << endl;

  try {
    ptime t = microsec_clock::universal_time();
    string datetime = to_iso_extended_string(t) + "Z";
    label += "_" + datetime;
  } catch (...) {
    cout << "WARNING: An error occurred while generating the date part of a label." << endl
         << "Although this shouldn't occur (something weird is going on with your system!)," << endl
         << "you can likely circumvent this issue by using the --output command line switch)" << endl
         << "to provide you own output label." << endl;
    if(n_rnd_char < 5)
      n_rnd_char = 5;
  }

  if (n_rnd_char > 0) {
    label += "_";
    for (size_t i = 0; i < n_rnd_char; i++)
      label += r_char(rng);
  }

  if (verbosity >= Verbosity::debug)
    cout << "Generated random label " << label << endl;

  return label;
}
