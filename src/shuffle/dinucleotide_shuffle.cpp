#include "dinucleotide_shuffle.hpp"
#include <cassert>
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <random>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

// based on altschulEriksonDinuclShuffle.py
// P. Clote, Oct 2003

typedef char Nucl;
typedef pair<Nucl,Nucl> Edge;
typedef list<Nucl> Nucls;
typedef list<Edge> EdgeList;
typedef map<Nucl,size_t> NuclCount;
typedef map<Edge,size_t> DinuclCount;

typedef map<Nucl,vector<Nucl>> NuclList;

const string nuclList   = "ACGTN";

void computeCountAndLists(const string &s, NuclCount &nuclCnt, DinuclCount &dinuclCnt, NuclList &nl) {
  // Compute count and lists
  nuclCnt[s[0]] = 1;
  size_t nuclTotal = 1;
  size_t dinuclTotal = 0;
  for(size_t i = 0; i < s.size() - 1; i++) {
    Nucl x = s[i]; Nucl y = s[i+1];
    nl[x].push_back(y);
    nuclCnt[y] += 1; nuclTotal  += 1;
    dinuclCnt[make_pair(x,y)] += 1; dinuclTotal += 1;
  }
  assert(nuclTotal==s.size());
  assert(dinuclTotal==s.size()-1);
}

template <class T> struct Eulerian {

  bool valid;
  EdgeList edges;
  Nucls nucls;
  char lastChar;
  std::uniform_real_distribution<double> r_unif;

  Eulerian() : valid(false), edges(), nucls(), lastChar('X') { };

  Eulerian(const string &s, T &rng) : valid(false), edges(), nucls(), lastChar(s.back()), r_unif(0,1) {
    NuclCount nuclCnt;
    DinuclCount dinuclCnt;
    NuclList nl;
    computeCountAndLists(s, nuclCnt, dinuclCnt, nl);

    // compute nucleotides appearing in s
    for(Nucl x: "ACGTN") if(s.find(x) != string::npos) nucls.push_back(x);

    for(Nucl x: nucls) if(x != lastChar) edges.push_back(make_pair(x, chooseEdge(x,dinuclCnt, rng)));

    valid = connectedToLast();
  };
  private:
  bool connectedToLast() const {
    map<Nucl,bool> D;
    for(auto x: nucls) D[x] = false;
    for(auto &edge: edges) if(edge.second == lastChar) D[edge.first] = true;
    for(size_t i = 0; i < 3; i++)
      for(auto &edge: edges) if(D[edge.second] == true) D[edge.first] = true;
    for(Nucl x: nucls) if(x != lastChar and D[x] == false) return false;
    return true;
  };

  Nucl chooseEdge(Nucl x, DinuclCount dinuclCnt, T &rng) {
    double z = r_unif(rng);
    double denom = 0;
    for(auto y: "ACGTN") {
      auto iter = dinuclCnt.find({x,y});
      if(iter != end(dinuclCnt))
        denom += iter->second;
    }
    double numerator = 0;
    for(auto y: "ACGT") {
      auto iter = dinuclCnt.find({x, y});
      if(iter != end(dinuclCnt))
        numerator += iter->second;
      if(z < numerator / denom) {
        dinuclCnt[{x,y}] -= 1;
        return y;
      }
    }
    dinuclCnt[{x,'N'}] -= 1;
    return 'N';
  };
};


string dinucleotideShuffle(const string &s_, size_t seed) {
  mt19937 rng;
  rng.seed(seed);
  string s(s_);
  boost::algorithm::to_upper(s);
  Eulerian<mt19937> eulerian;
  while(not eulerian.valid)
    eulerian = Eulerian<mt19937>(s, rng);

  NuclCount nuclCnt;
  DinuclCount dinuclCnt;
  NuclList nl;
  computeCountAndLists(s, nuclCnt, dinuclCnt, nl);

  // remove last edges from each vertex list, shuffle, then add back
  // the removed edges at end of vertex lists.
  for(Edge x: eulerian.edges) {
    auto &iter = nl.find(x.first)->second;
    iter.erase(std::find(begin(iter), end(iter), x.second));
  }
  for(Nucl x: eulerian.nucls) shuffle(begin(nl[x]), end(nl[x]), rng);
  for(Edge x: eulerian.edges) nl[x.first].push_back(x.second);

  // construct the Eulerian path
  string l = s.substr(0,1);
  Nucl prevCh = s[0];
  for(size_t i = 0; i < s.size() - 2; i++) {
    Nucl ch = nl[prevCh][0];
    l += ch;
    nl[prevCh].erase(begin(nl[prevCh]));
    prevCh = ch;
  }
  l += s.back();
  return l;
}

