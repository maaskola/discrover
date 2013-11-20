/*
 * =====================================================================================
 *
 *       Filename:  specification.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/22/2012 02:09:19 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */

#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include "../aux.hpp"
#include "specification.hpp"

using namespace std;

string readfile(const string &path) {
  string s;
  ifstream ifs(path.c_str());
  while(ifs.good()) {
    char c = ifs.get();
    if(ifs.good())
      s += c;
  }
  return(s);
}

namespace Specification {
  DataSet::DataSet() : series(""), path(""), is_shuffle(false), motifs() {
    if(false)
      cout << "Specification::DataSet standard constructor" << endl;
  };
  DataSet::DataSet(const DataSet &spec) : series(spec.series), path(spec.path), is_shuffle(spec.is_shuffle), motifs(spec.motifs) {
    if(false)
      cout << "Specification::DataSet copy constructor: '" << spec.path << "'" << endl;
  };
  DataSet::DataSet(const string &token, bool is_shuffle_) : DataSet() {
    if(false)
      cout << "Specification::DataSet string constructor: '" << token << "'." << endl;
    // the format is [NAMES:[SERIES:]]seq.fa e.g. pum,qki:bg0:seq.fa or :bg0:seq.fa or seq.fa
    // so if a file name contains at least one ":" then one would just have to write ::file_name_with_:_.fa
    is_shuffle = is_shuffle_;
    size_t pos;
    if((pos = token.find(":")) != string::npos) {
      string motif_token = token.substr(0, pos);
      string rest = token.substr(pos+1);
      for(auto &motif: tokenize(motif_token, ","))
        motifs.insert(motif);
      if((pos = rest.find(":")) != string::npos) {
        series = rest.substr(0, pos);
        path = rest.substr(pos+1);
      } else
        path = rest;
    } else
      path = token;

    if(boost::filesystem::exists(path)) {
      if(not boost::filesystem::is_regular_file(path)) {
        cout << "Error: FASTA file " << path << " is not a regular file." << endl;
        exit(-1);
      }
    } else {
      cout << "Error: FASTA file " << path << " does not exist." << endl;
      exit(-1);
    }

    if(false) {
      cout << "Constructed DataSet:\npath = " << path << "\n" << "series = " << series << "\n" << "motifs  =";
      for(auto x: motifs) cout << " " << x;
      cout << endl;
      cout << "This is" << (is_shuffle ? "" : " not") << " a shuffle." << endl;
    }
  }


  Motif::Motif(const Motif &s) : kind(s.kind), specification(s.specification), name(s.name), insertions(s.insertions), lengths(s.lengths)
  {
  }

  Motif::Motif(const string &s) : kind(Motif::Kind::Seed), specification(s), name(""), insertions(), lengths()
  {
    // cout << "Constructing Specification::Motif from '" << s << "'." << endl;
    size_t pos;
    if((pos = specification.find(":")) != string::npos) {
      name = specification.substr(0,pos);
      specification = specification.substr(pos+1);
    }
    if((pos = specification.find(":")) != string::npos) {
      insertions = parse_list(specification.substr(0,pos));
      specification = specification.substr(pos+1);
    } else {
      insertions = vector<size_t>();
    }
    if(boost::filesystem::exists(specification) and boost::filesystem::is_regular_file(specification))
      kind = Kind::File;
    else {
      if(specification.find_first_not_of("acgtubdhvkmwsrynACGTUBDHVKMWSRYN") == string::npos)
        kind = Kind::Seed;
      else {
        kind = Kind::Plasma;
        lengths = parse_list(specification);
      }
    }
    if(false) {
      cout << "kind = '" << static_cast<int>(kind) << "' name = '" << name << "' spec = '" << specification << "' ins =";
      for(auto x: insertions)
        cout << " " << x;
      cout << " lengths =";
      for(auto x: lengths)
        cout << " " << x;
      cout << endl;
    }
  }
}

namespace Specification {
  namespace Series {
    Expression make_series(const string &path) {
      // cout << "make_series(" << path << ")" << endl;
      Expression expr;
      string s = path;
      while(s.size() > 0) {
        size_t pos = s.find_last_of("+-");
        // cout << "pos = " << pos << endl;
        if(pos == string::npos) {
          // cout << "last one: " << s << endl;
          Atom atom = {+1, s};
          expr.push_back(atom);
          break;
        } else {
          Atom atom = {(s[pos] == '+' ? +1 : -1), s.substr(pos+1)};
          // cout << "next one: " << s.substr(pos+1) << endl;
          expr.push_back(atom);
          s = s.substr(0, pos);
        }
      }
      reverse(begin(expr), end(expr));
      return(expr);
    }
    string to_string(const Atom &atom) {
      return(string() + (atom.sign > 0 ? "+" : "-") + atom.series);
    }
    ostream &operator<<(ostream &out, const Atom &atom) {
      out << to_string(atom);
      return(out);
    }

    std::string to_string(const Expression &expr) {
      bool first = true;
      string s;
      for(auto &x: expr) {
        if(first)
          first = false;
        else
          s += ",";
        s += to_string(x);
      }
      return(s);
    }

    std::ostream &operator<<(std::ostream &out, const Expression &expr) {
      out << to_string(expr);
      return(out);
    }
  }
  istream &operator>>(istream &in, Specification::DataSet &spec) {
    // the format is [0,1,2:]seq.fa e.g. 0,1,2:seq.fa or 0,1:seq.fa or seq.fa
    // so if a file name contains at least one ":" then one would just have to write :file_name_with_:_.fa
    string token;
    in >> token;
    spec = Specification::DataSet(token);
    return(in);
  }

  istream &operator>>(istream &in, Specification::Motif &m) {
    string token;
    in >> token;
    m = Specification::Motif(token);
    return(in);
  }

  std::ostream &operator<<(std::ostream &out, const Motif &m) {
    out << to_string(m);
    return(out);
  }

  std::ostream &operator<<(std::ostream &out, const DataSet &spec) {
    out << to_string(spec);
    return(out);
  }

  string to_string(const Motif &spec) {
    string s = spec.name + ":";
    bool first = true;
    for(auto &l: spec.insertions) {
      if(first)
        first = false;
      else
        s += ",";
      s += boost::lexical_cast<string>(l);
    }
    s += spec.specification;
    return(s);
  }
  string to_string(const DataSet &spec) {
    string s = "";
    bool first = true;
    for(auto &m: spec.motifs) {
      if(first)
        first = false;
      else
        s += ",";
      s += m;
    }
    s += ":" + spec.series + ":" + spec.path;
    if(spec.is_shuffle)
      s += "_SHUFFLE";
    return(s);
  }
}

