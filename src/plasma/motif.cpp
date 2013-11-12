/*
 * =====================================================================================
 *
 *       Filename:  motif.cpp
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

#include <iostream>
#include "aux.hpp"
#include "code.hpp"
#include "options.hpp"

using namespace std;

namespace Seeding {
  double information_content(const string &motif) {
    double ic = 0;
    for(auto &c: motif)
      switch(c) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
          ic += 2;
          break;
        case 'w':
        case 's':
        case 'm':
        case 'k':
        case 'r':
        case 'y':
          ic += 1;
          break;
        case 'b':
        case 'd':
        case 'h':
        case 'v':
          ic += log(4.0/3.0) / log(2.0);
          break;
        case 'n':
          ic += 0;
          break;
        default:
          break;
      }
    return(ic);
  }
  size_t motif_degeneracy(const string &motif) {
    size_t s = 0;
    for(auto &c: motif)
      switch(c) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
          s += 0;
          break;
        case 'w':
        case 's':
        case 'm':
        case 'k':
        case 'r':
        case 'y':
          s += 1;
          break;
        case 'b':
        case 'd':
        case 'h':
        case 'v':
          s += 2;
          break;
        case 'n':
          s += 3;
          break;
        default:
          break;
      }
    return(s);
  }

  list<string> all_generalizations(const string &motif) {
    list<string> generalizations;
    // cout << "Generalizations of " << to_string() << endl;
    string generalization = motif;
    for(size_t i = 0; i < motif.size(); i++) {
      char current = generalization[i];
      switch(current) {
        case 'a':
          generalization[i] = 'm';
          generalizations.push_back(generalization);
          generalization[i] = 'r';
          generalizations.push_back(generalization);
          generalization[i] = 'w';
          generalizations.push_back(generalization);
          break;
        case 'c':
          generalization[i] = 'm';
          generalizations.push_back(generalization);
          generalization[i] = 's';
          generalizations.push_back(generalization);
          generalization[i] = 'y';
          generalizations.push_back(generalization);
          break;
        case 'g':
          generalization[i] = 'r';
          generalizations.push_back(generalization);
          generalization[i] = 's';
          generalizations.push_back(generalization);
          generalization[i] = 'k';
          generalizations.push_back(generalization);
          break;
        case 't':
          generalization[i] = 'w';
          generalizations.push_back(generalization);
          generalization[i] = 'y';
          generalizations.push_back(generalization);
          generalization[i] = 'k';
          generalizations.push_back(generalization);
          break;
        case 'w':
          generalization[i] = 'd';
          generalizations.push_back(generalization);
          generalization[i] = 'h';
          generalizations.push_back(generalization);
        case 's':
          generalization[i] = 'b';
          generalizations.push_back(generalization);
          generalization[i] = 'v';
          generalizations.push_back(generalization);
        case 'm':
          generalization[i] = 'h';
          generalizations.push_back(generalization);
          generalization[i] = 'v';
          generalizations.push_back(generalization);
        case 'k':
          generalization[i] = 'b';
          generalizations.push_back(generalization);
          generalization[i] = 'd';
          generalizations.push_back(generalization);
        case 'r':
          generalization[i] = 'd';
          generalizations.push_back(generalization);
          generalization[i] = 'v';
          generalizations.push_back(generalization);
        case 'y':
          generalization[i] = 'b';
          generalizations.push_back(generalization);
          generalization[i] = 'h';
          generalizations.push_back(generalization);
        case 'b':
        case 'd':
        case 'h':
        case 'v':
          generalization[i] = 'n';
          generalizations.push_back(generalization);
        default:
          break;
      }
      generalization[i] = current;
    }
    return(generalizations);
  }
}


