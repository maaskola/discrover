/*
 * =====================================================================================
 *
 *       Filename:  code.cpp
 *
 *    Description:  
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <limits>
#include "code.hpp"

using namespace std;

namespace Seeding {
  vector<symbol_t> construct_code() {
    const size_t table_size = numeric_limits<symbol_t>::max() + 1;
    vector<symbol_t> table(table_size, 0);
    for(symbol_t i = 0; i < Symbol.size(); i++)
      table[static_cast<size_t>(Symbol[i])] = i;
    for(symbol_t x = 'A'; x <= 'Z'; x++)
      table[static_cast<size_t>(x)] = table['a' + x - 'A'];
    //  for(symbol_t i = 0; i < 128; i++)
    //    cout << i << " " << static_cast<size_t>(table[i]) << endl;
    return(table);
  };

  const vector<symbol_t> Code = construct_code();

  string iupac2regex(const string &s) {
    string r;
    for(auto c: s)
      switch(tolower(c)) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'u':
          r += c;
          break;
        case 'k':
          r += "[gt]";
          break;
        case 'm':
          r += "[ac]";
          break;
        case 's':
          r += "[cg]";
          break;
        case 'w':
          r += "[at]";
          break;
        case 'r':
          r += "[ag]";
          break;
        case 'y':
          r += "[ct]";
          break;
        case 'd':
          r += "[agt]";
          break;
        case 'b':
          r += "[cgt]";
          break;
        case 'h':
          r += "[act]";
          break;
        case 'v':
          r += "[acg]";
          break;
        case 'n':
          r += ".";
          break;
      }
    return(r);
  }

  bool iupac_included(char r, char q)
  {
    q = tolower(q);
    r = tolower(r);
    if(q==r)
      return(true);
    switch(r) {
      case 'a':
        switch(q) {
          case 'w': case 'm': case 'r': case 'd': case 'h': case 'v': case 'n':
            return true;
          default:
            return false;
        }
      case 'c':
        switch(q) {
          case 's': case 'm': case 'y': case 'b': case 'h': case 'v': case 'n':
            return true;
          default:
            return false;
        }
      case 'g':
        switch(q) {
          case 's': case 'k': case 'r': case 'b': case 'd': case 'v': case 'n':
            return true;
          default:
            return false;
        }
      case 't':
        switch(q) {
          case 'u': case 'w': case 'k': case 'y': case 'b': case 'd': case 'h': case 'n':
            return true;
          default:
            return false;
        }
      case 'u':
        switch(q) {
          case 't': case 'w': case 'k': case 'y': case 'b': case 'd': case 'h': case 'n':
            return true;
          default:
            return false;
        }
      default:
        return(false);
    }
  }
}

seq_type encode(const string &s) {
  size_t n = s.size();
  seq_type vec(n);
  auto iter = begin(s);
  for(auto &v: vec)
    v = Seeding::Code[*iter++];
  return vec;
}

std::string decode(const seq_type &seq) {
  size_t n = seq.size();
  string s(n, ' ');
  auto iter = begin(seq);
  for(auto &c: s)
    c = Seeding::Symbol[*iter++];
  return s;
}

void add_sequence(vector<symbol_t> &s, const string &seq) {
  for(auto x: seq)
    s.push_back(Seeding::Code[x]);
}

using nucl_vector_type = vector<bool>;

nucl_vector_type build_pure_nucl_vector() {
  const size_t n = 16;
  nucl_vector_type v(n);
  for (size_t i = 0; i < n; ++i) v[i] = false;
  for (auto nucl : encode("acgt")) v[nucl] = true;
  return v;
}

const static nucl_vector_type pure_nucl_vector = build_pure_nucl_vector();

bool pure_nucleotide(symbol_t s) {
  return pure_nucl_vector[s];
}

bool degenerate_nucleotide(symbol_t s) {
  return not pure_nucl_vector[s];
}
