/*
 * =====================================================================================
 *
 *       Filename:  motif.cpp
 *
 *    Description:  
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include "../aux.hpp"
#include "code.hpp"
#include "motif.hpp"
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

  vector<vector<seq_type::value_type>> build_generalization_table() {
    using val_t = seq_type::value_type;
    vector<vector<val_t>> t;
    for(val_t s = 0; s < 16; ++s) {
      vector<val_t> v;
      for(size_t i = 0; i < 4; ++i) {
        val_t x = s | (1 << i);
        if(x != s)
          v.push_back(x);
      }
      t.push_back(v);
    }
    if(false) {
      cout << "generalization table:" << endl;
      for(val_t i = 0; i < t.size(); ++i) {
        cout << Seeding::Symbol[i] << ":";
        for(auto &x: t[i])
          cout << " " << int(x) << ":"
            << "'" << Seeding::Symbol[x] << "'";
        cout << endl;
      }
    }
    return t;
  }

  const vector<vector<seq_type::value_type>> generalization_table = build_generalization_table();

  vector<seq_type> all_generalizations(const seq_type &motif) {
    vector<seq_type> generalizations;
    seq_type generalization = motif;
    // cout << "Generalizations of " << decode(generalization) << endl;
    for(size_t i = 0; i < motif.size(); i++) {
      auto current = generalization[i];
      for(auto &x: generalization_table[current]) {
        generalization[i] = x;
        // cout << "Generalization -> " << decode(generalization) << endl;
        generalizations.push_back(generalization);
      }
      generalization[i] = current;
    }
    return(generalizations);
  }
}


