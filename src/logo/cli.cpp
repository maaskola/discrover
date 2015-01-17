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
#include "../aux.hpp"

using namespace std;

namespace Logo {
istream &operator>>(istream &is, Type &type) {
  string Word;
  is >> Word;
  string word = string_tolower(Word);
  if (word == "seq" or word == "sequence" or word == "info"
      or word == "information")
    type = Type::Sequence;
  else if (word == "freq" or word == "frequency")
    type = Type::Frequency;
  else
    throw Exception::InvalidType(Word);
  return is;
}
istream &operator>>(istream &is, Alphabet &alphabet) {
  string Word;
  is >> Word;
  string word = string_tolower(Word);
  if (word == "rna")
    alphabet = Alphabet::RNA;
  else if (word == "dna")
    alphabet = Alphabet::DNA;
  else if (word == "undef" or word == "undefined")
    alphabet = Alphabet::Undefined;
  else
    throw Exception::InvalidAlphabet(Word);

  return is;
}
istream &operator>>(istream &is, Order &order) {
  string Word;
  is >> Word;
  string word = string_tolower(Word);
  if (word == "alpha" or word == "alphabetic")
    order = Order::Alphabetic;
  else if (word == "freq" or word == "frequency")
    order = Order::Frequency;
  else
    throw Exception::InvalidOrder(Word);
  return is;
}
istream &operator>>(istream &is, Palette &palette) {
  string Word;
  is >> Word;
  string word = string_tolower(Word);
  if (word == "default")
    palette = Palette::Default;
  else if (word == "solarized")
    palette = Palette::Solarized;
  else if (word == "tetrad")
    palette = Palette::Tetrad;
  else
    throw Exception::InvalidPalette(Word);
  return is;
}

namespace Exception {
InvalidType::InvalidType(const string &token)
    : runtime_error("Error: invalid logo type '" + token + "'.\n"
                    + "Available are: 'seq' and 'freq'.") {}
InvalidAlphabet::InvalidAlphabet(const string &token)
    : runtime_error( "Error: invalid alphabet '" + token + "'.\n"
    "Available are: 'DNA' and 'RNA'."){
    }
InvalidOrder::InvalidOrder(const string &token)
    : runtime_error("Error: invalid order type '" + token + "'.\n"
                    + "Available are: 'alpha' and 'freq'.") {}
InvalidPalette::InvalidPalette(const string &token)
    : runtime_error("Error: invalid color palette '" + token + "'.\n"
                    + "Available are: 'default', 'solarized', and 'tetrad'.") {}
}  // namespace Exception
}  // namespace Logo

boost::program_options::options_description gen_logo_options_description(
    Logo::Options &options, Logo::CLI mode, size_t cols, const string &name) {
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
    ("logo", po::value(&options.type)->default_value(Logo::Type::Sequence, "info"), "Which kind of logo to create; 'info' for information-type sequence logo (position height scaled by information content), 'freq' for frequency logo.")
    ("alphabet", po::value(&options.alphabet), "Which alphabet to use; can be either 'RNA' or 'DNA'. If left unspecified, 'DNA' is chosen if --revcomp is used, and 'RNA' otherwise.")
    ("order", po::value(&options.order)->default_value(Logo::Order::Frequency, "freq"), "How to vertically order the nucleotides; can be either 'alpha' for alphabetic order or 'freq' for most frequent at top.")
    ("pal", po::value(&options.palette)->default_value(Logo::Palette::Default, "default"), "Color palette to use; available are 'default', 'solarized', 'tetrad'.")
    ("scale", po::value(&options.scale)->default_value(100.0, "100"), "Height in pixels of the nucleotide stacks in the sequence logos.")
    ;
  if (mode != Logo::CLI::HMM)
    desc.add_options()
      ("absent", po::value(&options.absent)->default_value(0.03, "0.03"), "Use this frequency for absent nucleotides when creating logos for IUPAC regular expression motifs.")
      ;
  else
    options.absent = 0.03;

  if (mode == Logo::CLI::Full)
    desc.add_options()
      ("revcomp,r", po::bool_switch(&options.revcomp), "Generate sequence logos for forward and reverse complementary strand.")
        ;
  else
    options.revcomp = false;

  return desc;
}
