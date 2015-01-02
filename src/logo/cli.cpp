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
#include "options.hpp"
#include "../executioninformation.hpp"

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
  alphabet = Alphabet::RNA;
  return is;
}
istream &operator>>(istream &is, Order &order) {
  order = Order::Frequency;
  return is;
}
istream &operator>>(istream &is, Palette &palette) {
  palette = Palette::Default;
  return is;
}
}

boost::program_options::options_description gen_logo_options_description(
    Logo::Options &options, bool iupac_mode, size_t cols,
    const string &name) {
  namespace po = boost::program_options;
  po::options_description desc(name, cols);
  if (iupac_mode)
    desc.add_options()
      ("pdf", po::bool_switch(&options.pdf_logo), "Generate PDF files with sequence logos of the found motifs.")
      ("png", po::bool_switch(&options.png_logo), "Generate PDF files with sequence logos of the found motifs.")
      ;
  else
    desc.add_options()
      ("nopdf", po::bool_switch()->notifier([&](bool val){options.pdf_logo = not val;}), "Generate PDF files with sequence logos of the found motifs.")
      ("nopng", po::bool_switch()->notifier([&](bool val){options.png_logo = not val;}), "Generate PNG files with sequence logos of the found motifs.")
      ;
  desc.add_options()
    ("axes", po::bool_switch(&options.axes), "Include axes in sequence logos.")
    ("logotype", po::value<Logo::Type>(&options.type)->default_value(Logo::Type::Sequence, "seq"), "Which kind of logo to create; 'seq' for sequence logo (position height scaled by information content), 'freq' for frequency logo.")
    ;
  if (iupac_mode)
    desc.add_options()
      ("absent", po::value<double>(&options.absent)->default_value(0.03, "0.03"), "Use this frequency for absent nucleotides when creating logos for IUPAC regular expression motifs.")
      ;
  else
    options.absent = 0.03;

  return (desc);
}
