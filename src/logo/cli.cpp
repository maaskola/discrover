/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Created:  Fri Jan 02 20:32:15 2014 +0100
 *         Author:  Jonas Maaskola (jonas@maaskola.de)
 *
 * =====================================================================================
 */

#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include <iostream>
#include "../executioninformation.hpp"
#include "cli.hpp"

using namespace std;

namespace Logo {
istream &operator>>(istream &is, Type &type) {
  string word;
  is >> word;
  if (word == "seq" or word == "sequence")
    type = Type::Sequence;
  else if (word == "freq" or word == "frequency")
    type = Type::Frequency;
  else {
    cout << "Can't parse logo type " << word << "." << endl;
    exit(-1);
  }
  return is;
}
istream &operator>>(istream &is, Alphabet &alphabet) {
  string word;
  is >> word;
  if (word == "rna" or word == "RNA")
    alphabet = Alphabet::RNA;
  else if (word == "dna" or word == "DNA")
    alphabet = Alphabet::DNA;
  else if (word == "undef" or word == "undefined" or word == "Undef" or word == "Undefined")
    alphabet = Alphabet::Undefined;
  else {
    cout << "Can't parse alphabet " << word << "." << endl;
    exit(-1);
  }

  return is;
}
istream &operator>>(istream &is, Order &order) {
  string word;
  is >> word;
  if (word == "alpha" or word == "alphabetic")
    order = Order::Alphabetic;
  else if (word == "freq" or word == "frequency")
    order = Order::Frequency;
  else {
    cout << "Can't parse order type " << word << "." << endl;
    exit(-1);
  }
  return is;
}
istream &operator>>(istream &is, Palette &palette) {
  palette = Palette::Default;
  return is;
}
}

boost::program_options::options_description gen_logo_options_description(
    Logo::Options &options, Logo::CLI mode, size_t cols,
    const string &name) {
  namespace po = boost::program_options;
  po::options_description desc(name, cols);
  if (mode == Logo::CLI::IUPAC)
    desc.add_options()
      ("pdf", po::bool_switch(&options.pdf_logo), "Generate PDF files with sequence logos of the found motifs.")
      ("png", po::bool_switch(&options.png_logo), "Generate PDF files with sequence logos of the found motifs.")
      ;
  else
    desc.add_options()
      ("nopdf", po::bool_switch()->notifier([&](bool val){options.pdf_logo = not val;}), "Do not generate PDF files with sequence logos of the found motifs.")
      ("nopng", po::bool_switch()->notifier([&](bool val){options.png_logo = not val;}), "Do not generate PNG files with sequence logos of the found motifs.")
      ;
  desc.add_options()
    ("axes", po::bool_switch(&options.axes), "Include axes in sequence logos.")
    ("logotype", po::value<Logo::Type>(&options.type)->default_value(Logo::Type::Sequence, "seq"), "Which kind of logo to create; 'seq' for sequence logo (position height scaled by information content), 'freq' for frequency logo.")
    ("alphabet", po::value<Logo::Alphabet>(&options.alphabet), "Which alphabet to use; can be either 'RNA' or 'DNA'. If left unspecified, then 'DNA' is chosen if --revcomp is used, and otherwise 'RNA'.")
    ("order", po::value<Logo::Order>(&options.order)->default_value(Logo::Order::Frequency, "freq"), "How to vertically order the nucleotides; can be either 'alpha' for alphabetic order or 'freq' for most frequent at top.")
    ;
  if (mode == Logo::CLI::IUPAC)
    desc.add_options()
      ("absent", po::value<double>(&options.absent)->default_value(0.03, "0.03"), "Use this frequency for absent nucleotides when creating logos for IUPAC regular expression motifs.")
      ;
  else
    options.absent = 0.03;

  if (mode == Logo::CLI::Full)
    desc.add_options()
      ("rc", po::bool_switch(&options.revcomp), "Generate sequence logos for forward and reverse complementary strand.")
        ;
  else
    options.revcomp = false;

  return (desc);
}
