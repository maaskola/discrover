/*
 * =====================================================================================
 *
 *       Filename:  code.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31.05.2012 06:47:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *   Organization:  
 *
 * =====================================================================================
 */

#include "code.hpp"

using namespace std;

namespace Seeding {

  char* construct_code() {
    char* table = new char[128];
    for(size_t i = 0; i < 128; i++)
      table[i] = 0;
    for(size_t i = 0; i < 16; i++)
      table[Symbol[i]] = i;
    table['A'] = table['a'];
    table['C'] = table['c'];
    table['G'] = table['g'];
    table['T'] = table['t'];
    table['U'] = table['u'];
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

}

