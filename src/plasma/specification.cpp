/*
 * =====================================================================================
 *
 *       Filename:  specification.cpp
 *
 *    Description:
 *
 *        Created:  06/22/2012 02:09:19 PM
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include <boost/filesystem.hpp>
#include "../aux.hpp"
#include "specification.hpp"
#include "io.hpp"

using namespace std;

namespace Specification {

Set::Set() : contrast(""), path(""), is_shuffle(false), motifs() {
  if (false)
    cout << "Specification::Set standard constructor" << endl;
};

Set::Set(const Set &spec)
    : contrast(spec.contrast),
      path(spec.path),
      is_shuffle(spec.is_shuffle),
      motifs(spec.motifs) {
  if (false)
    cout << "Specification::Set copy constructor: '" << spec.path << "'"
         << endl;
};

Set::Set(const string &token, bool is_shuffle_) : Set() {
  if (false)
    cout << "Specification::Set string constructor: '" << token << "'." << endl;
  // the format is [NAMES:[SERIES:]]seq.fa e.g. pum,qki:bg0:seq.fa or
  // :bg0:seq.fa or seq.fa
  // so if a file name contains at least one ":" then one would just have to
  // write ::file_name_with_:_.fa
  is_shuffle = is_shuffle_;
  size_t pos;
  if ((pos = token.find(":")) != string::npos) {
    string motif_token = token.substr(0, pos);
    string rest = token.substr(pos + 1);
    for (auto &motif : tokenize(motif_token, ","))
      motifs.insert(motif);
    if ((pos = rest.find(":")) != string::npos) {
      contrast = rest.substr(0, pos);
      path = rest.substr(pos + 1);
    } else
      path = rest;
  } else
    path = token;

  if (not boost::filesystem::exists(path))
    throw ::Exception::File::Existence(path);
  else if (not boost::filesystem::is_regular_file(path))
    throw ::Exception::File::NoRegularFile(path);

  if (false) {
    cout << "Constructed Set:\npath = " << path << "\n"
         << "contrast = " << contrast << "\n"
         << "motifs  =";
    for (auto x : motifs)
      cout << " " << x;
    cout << endl;
    cout << "This is" << (is_shuffle ? "" : " not") << " a shuffle." << endl;
  }
}

string Set::name() const {
  return (is_shuffle ? "shuffle(" : "") + path + (is_shuffle ? ")" : "");
}

Motif::Motif(const Motif &s)
    : kind(s.kind),
      specification(s.specification),
      name(s.name),
      insertions(s.insertions),
      lengths(s.lengths),
      multiplicity(s.multiplicity){}

Motif::Motif(const string &s)
    : kind(Motif::Kind::Seed),
      specification(s),
      name(""),
      insertions(),
      lengths(),
      multiplicity(1) {
  // cout << "Constructing Specification::Motif from '" << s << "'." << endl;
  size_t pos;
  if ((pos = specification.find(":")) != string::npos) {
    name = specification.substr(0, pos);
    specification = specification.substr(pos + 1);
  }
  if ((pos = specification.find(":")) != string::npos) {
    insertions = parse_list(specification.substr(0, pos));
    specification = specification.substr(pos + 1);
  } else {
    insertions = vector<size_t>();
  }
  if (boost::filesystem::exists(specification)
      and boost::filesystem::is_regular_file(specification))
    kind = Kind::File;
  else {
    if (specification.find_first_not_of("acgtubdhvkmwsrynACGTUBDHVKMWSRYN")
        == string::npos)
      kind = Kind::Seed;
    else {
      kind = Kind::Plasma;
      if ((pos = specification.find("x")) != string::npos) {
        multiplicity = atoi(specification.substr(pos + 1).c_str());
        specification = specification.substr(0, pos);
      }
      lengths = parse_list(specification);
    }
  }
  if (false) {
    cout << "kind = '" << static_cast<int>(kind) << "' name = '" << name
         << "' spec = '" << specification << "' ins =";
    for (auto x : insertions)
      cout << " " << x;
    cout << " lengths =";
    for (auto x : lengths)
      cout << " " << x;
    cout << endl;
    cout << "multiplicity = " << multiplicity << endl;
  }
}

namespace Contrast {

Expression make_contrast(const string &path) {
  // cout << "make_contrast(" << path << ")" << endl;
  Expression expr;
  string s = path;
  while (s.size() > 0) {
    size_t pos = s.find_last_of("+-");
    // cout << "pos = " << pos << endl;
    if (pos == string::npos) {
      // cout << "last one: " << s << endl;
      Atom atom = {+1, s};
      expr.push_back(atom);
      break;
    } else {
      Atom atom = {(s[pos] == '+' ? +1 : -1), s.substr(pos + 1)};
      // cout << "next one: " << s.substr(pos+1) << endl;
      expr.push_back(atom);
      s = s.substr(0, pos);
    }
  }
  reverse(begin(expr), end(expr));
  return expr;
}

string to_string(const Atom &atom) {
  return string() + (atom.sign > 0 ? "+" : "-") + atom.contrast;
}

ostream &operator<<(ostream &out, const Atom &atom) {
  out << to_string(atom);
  return out;
}

std::string to_string(const Expression &expr) {
  bool first = true;
  string s;
  for (auto &x : expr) {
    if (first)
      first = false;
    else
      s += ",";
    s += to_string(x);
  }
  return s;
}

std::ostream &operator<<(std::ostream &out, const Expression &expr) {
  out << to_string(expr);
  return out;
}
}

istream &operator>>(istream &in, Specification::Set &spec) {
  // the format is [0,1,2:]seq.fa e.g. 0,1,2:seq.fa or 0,1:seq.fa or seq.fa
  // so if a file name contains at least one ":" then one would just have to
  // write :file_name_with_:_.fa
  string token;
  in >> token;
  spec = Specification::Set(token);
  return in;
}

istream &operator>>(istream &in, Specification::Motif &m) {
  string token;
  in >> token;
  m = Specification::Motif(token);
  return in;
}

std::ostream &operator<<(std::ostream &out, const Motif &m) {
  out << to_string(m);
  return out;
}

std::ostream &operator<<(std::ostream &out, const Set &spec) {
  out << to_string(spec);
  return out;
}

string to_string(const Motif &spec) {
  string s = spec.name + ":";
  bool first = true;
  for (auto &l : spec.insertions) {
    if (first)
      first = false;
    else
      s += ",";
    s += std::to_string(l);
  }
  s += spec.specification;
  if (spec.multiplicity > 1)
    s += "x" + std::to_string(spec.multiplicity);
  return s;
}

string to_string(const Set &spec) {
  string s = "";
  bool first = true;
  for (auto &m : spec.motifs) {
    if (first)
      first = false;
    else
      s += ",";
    s += m;
  }
  s += ":" + spec.contrast + ":" + spec.path;
  if (spec.is_shuffle)
    s += "_SHUFFLE";
  return s;
}
}
