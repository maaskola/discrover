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

#include "code.hpp"

using namespace std;

namespace Seeding {
  char* Code = construct_code();

  char* construct_code() {
    char* table = new char[128];
    for(size_t i = 0; i < 128; i++)
      table[i] = 0;
    for(size_t i = 0; i < 16; i++)
      table[Symbol[i]] = i;
    for(char x = 'A'; x <= 'Z'; x++)
      table[x] = table['a' + x - 'A'];
    //  for(size_t i = 0; i < 128; i++)
    //    cout << i << " " << static_cast<int>(table[i]) << endl;
    return(table);
  };

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
