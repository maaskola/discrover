/*
 * =====================================================================================
 *
 *       Filename:  data.cpp
 *
 *    Description:
 *
 *        Created:  31.05.2012 06:47:48
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iostream>
#include <map>
#include <cstring>
#include "../aux.hpp"
#include "data.hpp"
#include "../verbosity.hpp"
#include "../sha1.hpp"

using namespace std;

string sha1hash(const string &s) {
  const char *t = s.c_str();
  unsigned char hash[20];
  char hexstring[41];  // 40 chars + a zero
  int end = (int)strlen(t);
  sha1::calc(t, end, hash);
  sha1::toHexString(hash, hexstring);
  string res;
  for (size_t i = 0; i < 40; i++)
    res += hexstring[i];
  return res;
}

namespace Data {
RemovalReport::RemovalReport(size_t n, size_t s)
    : nucleotides(0), sequences(0){};

RemovalReport &operator+=(RemovalReport &a, const RemovalReport &b) {
  a = a + b;
  return a;
}

RemovalReport operator+(const RemovalReport &a, const RemovalReport &b) {
  return RemovalReport(a.nucleotides + b.nucleotides,
                       a.sequences + b.sequences);
}
}

namespace Seeding {
Set::Set() : Data::Basic::Set<Fasta::Entry>() {}

Set::Set(const Specification::Set &s, bool revcomp, size_t n_seq)
    : Data::Basic::Set<Fasta::Entry>(s, revcomp, n_seq){};

Contrast::Contrast() : Data::Basic::Contrast<Set>() {}

Contrast::Contrast(const string &name, const Specification::Sets &paths,
                   bool revcomp, size_t n_seq)
    : Data::Basic::Contrast<Set>(name, paths, revcomp, n_seq) {}

Collection::Collection(const Specification::Sets &paths, bool revcomp,
                       size_t n_seq)
    : Data::Basic::Collection<Contrast>(paths, revcomp, n_seq) {}
}
