#include <iostream>
#include <fstream>
#include <cstdlib>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../terminal.hpp"
#include "../verbosity.hpp"
#include "../executioninformation.hpp"
#include "../aux.hpp"
#include "../plasma/code.hpp"
#include "../hmm/hmm.hpp"
#include "logo.hpp"
#include "cli.hpp"

using namespace std;

string gen_usage_string(const string &name) {
  string usage =
    "Generate sequence logos.\n"
    "\n"
    "Usage:\n"
    "  " + name + " path.hmm ... [options]\n"
    "\n"
    "Note that multiple paths can be given.";
  return usage;
}

Logo::matrix_t read_matrix(const string &path) {
  Logo::matrix_t matrix;
  ifstream ifs(path);
  string line;
  while (getline(ifs, line)) {
    Logo::column_t col(4, 0);
    stringstream ss(line);
    for (size_t i = 0; i < 4; i++)
      ss >> col[i];
    matrix.push_back(col);
  }
  return matrix;
}

Logo::matrix_t build_matrix(const string &motif, double absent) {
  const string nucls = "acgt";
  Logo::matrix_t matrix;
  for (auto pos : motif) {
    Logo::column_t col(4, 0);
    double z = 0;
    for (size_t i = 0; i < nucls.size(); i++)
      if (Seeding::iupac_included(nucls[i], pos)) {
        col[i] = 1;
        z++;
      }
    if (z == 0)
      for (size_t i = 0; i < nucls.size(); i++)
        col[i] = 1.0 / 4;
    else {
      for (size_t i = 0; i < nucls.size(); i++)
        if (col[i] > 0)
          col[i] = (1.0 - (4 - z) * absent) / z;
        else
          col[i] = absent;
    }
    matrix.push_back(col);
  }
  return matrix;
}

void draw_logos(const HMM &hmm, const Logo::Options &options,
                const string &label, size_t &motif_idx) {
  for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
    if (hmm.is_motif_group(group_idx)) {
      const string nucls = "acgt";
      Logo::matrix_t matrix;
      for (auto state : hmm.groups[group_idx].states) {
        Logo::column_t col(4, 0);
        for (size_t i = 0; i < nucls.size(); i++)
          col[i] = hmm.emission(state, i);
        matrix.push_back(col);
      }
      Logo::draw_logo(matrix, label + ".motif" + to_string(motif_idx++),
                      options);
      motif_idx++;
    }
}

int main(int argc, const char **argv) {
  const string default_error_msg
      = "Please inspect the command line help with -h or --help.";

  namespace po = boost::program_options;

  Logo::Options options;

  static const size_t MIN_COLS = 60;
  static const size_t MAX_COLS = 80;
  size_t cols = get_terminal_width();
  if (cols < MIN_COLS)
    cols = MIN_COLS;
  if (cols > MAX_COLS)
    cols = MAX_COLS;

  std::vector<string> hmm_paths, matrix_paths, iupacs;
  string label;
  ExecutionInformation exec_info
      = generate_exec_info(argv[0], GIT_DESCRIPTION, cmdline(argc, argv));

  po::options_description basic_options("Basic options", cols);

  basic_options.add_options()
    ("hmm", po::value<vector<string>>(&hmm_paths), "Path to .hmm file. May be given multiple times. Note: free arguments are also interpreted as .hmm files.")
    ("iupac,i", po::value<vector<string>>(&iupacs), "A motif given as a IUPAC regular expression. May be given multiple times.")
    ("matrix,m", po::value<vector<string>>(&matrix_paths), "Path to a file with a motif in matrix form. May be given multiple times.")
    ("output,o", po::value<string>(&label), string("Output file names are generated from this label. If this option is not specified the output label will be '" + exec_info.program_name + "_XXX' where XXX is a string to make the label unique.").c_str())
    ("help,h", "Produce help message.")
    ("version", "Print out the version.")
    ;

  po::options_description logo_options
      = gen_logo_options_description(options, Logo::CLI::Full, cols);

  po::options_description all_options;
  all_options.add(basic_options).add(logo_options);

  po::positional_options_description pos;
  pos.add("hmm", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(pos)
                  .run(),
              vm);
  } catch (po::unknown_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " not known." << endl << default_error_msg
         << endl;
    return EXIT_FAILURE;
  } catch (po::ambiguous_option &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " is ambiguous." << endl << default_error_msg
         << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_values &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::multiple_occurrences &e) {
    cout << "Error while parsing command line options:" << endl << "Option "
         << e.get_option_name() << " was specified multiple times." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_option_value &e) {
    cout << "Error while parsing command line options:" << endl
         << "The value specified for option " << e.get_option_name()
         << " has an invalid format." << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::too_many_positional_options_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Too many positional options were specified." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_syntax &e) {
    cout << "Error while parsing command line options:" << endl
         << "Invalid command line syntax." << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::invalid_command_line_style &e) {
    cout << "Error while parsing command line options:" << endl
         << "There is a programming error related to command line style."
         << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::reading_file &e) {
    cout << "Error while parsing command line options:" << endl
         << "The config file can not be read." << endl << default_error_msg
         << endl;
    return EXIT_FAILURE;
  } catch (po::validation_error &e) {
    cout << "Error while parsing command line options:" << endl
         << "Validation of option " << e.get_option_name() << " failed." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  } catch (po::error &e) {
    cout << "Error while parsing command line options:" << endl
         << "No further information as to the nature of this error is "
            "available, please check your command line arguments." << endl
         << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  if (vm.count("version") and not vm.count("help")) {
    cout << exec_info.program_name << " " << exec_info.hmm_version
         << " [" << GIT_BRANCH << " branch]" << endl;
    cout << GIT_SHA1 << endl;
    return EXIT_SUCCESS;
  }

  if (vm.count("help")) {
    cout << exec_info.program_name << " " << exec_info.hmm_version << endl;
    cout << "Copyright (C) 2015 Jonas Maaskola\n"
            "Provided under GNU General Public License Version 3 or later.\n"
            "See the file COPYING provided with this software for details of "
            "the license.\n" << endl;
    cout << limit_line_length(gen_usage_string(exec_info.program_name), cols)
         << endl << endl;
    cout << all_options << endl;
    return EXIT_SUCCESS;
  }

  try {
    po::notify(vm);
  } catch (po::required_option &e) {
    cout << "Error while parsing command line options:" << endl
         << "The required option " << e.get_option_name()
         << " was not specified." << endl << default_error_msg << endl;
    return EXIT_FAILURE;
  }

  // generate an output path stem if the user did not specify one
  if (not vm.count("output")) {
    label = generate_random_label(exec_info.program_name, 0, Verbosity::info);
    while (boost::filesystem::exists(label + ".hmm"))
      label = generate_random_label(exec_info.program_name, 5, Verbosity::info);
    cout << "Using \"" << label << "\" as label to generate output file names."
         << endl;
  }

  size_t motif_idx = 0;
  for (auto iupac : iupacs) {
    cout << "=> IUPAC motif " << iupac << endl;
    Logo::matrix_t matrix = build_matrix(iupac, options.absent);
    Logo::draw_logo(matrix, label + ".motif" + to_string(motif_idx++), options);
  }

  for (auto path : matrix_paths) {
    cout << "=> matrix file " << path << endl;
    Logo::matrix_t matrix = read_matrix(path);
    Logo::draw_logo(matrix, label + ".motif" + to_string(motif_idx++), options);
  }

  for (auto path : hmm_paths) {
    cout << "=> .hmm file " << path << endl;
    try {
      HMM hmm(path, Verbosity::info);
      draw_logos(hmm, options, label, motif_idx);
    } catch (runtime_error &e) {
      cout << e.what() << endl;
      return EXIT_FAILURE;
    }
  }

  if (motif_idx == 0)
    cout << "Nothing to do! " << default_error_msg << endl;

  return EXIT_SUCCESS;
}
