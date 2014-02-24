/*
 * =====================================================================================
 *
 *       Filename:  data.cpp
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
#include <map>
#include "../aux.hpp"
#include "data.hpp"
#include "../verbosity.hpp"
#include "../sha1.hpp"

using namespace std;

string sha1hash(const string &s)
{
  const char *t = s.c_str();
  unsigned char hash[20];
  char hexstring[41]; // 40 chars + a zero
  int end = (int) strlen(t);
  sha1::calc(t, end, hash);
  sha1::toHexString(hash, hexstring);
  string res;
  for(size_t i = 0; i < 40; i++)
    res += hexstring[i];
  return(res);
}

namespace Data {
  RemovalReport::RemovalReport(size_t n, size_t s) : nucleotides(0), sequences(0) { };
  RemovalReport &operator+=(RemovalReport &a, const RemovalReport &b) {
    a = a + b;
    return(a);
  }
  RemovalReport operator+(const RemovalReport &a, const RemovalReport &b) {
    return(RemovalReport(a.nucleotides + b.nucleotides, a.sequences + b.sequences));
  }
}


namespace Seeding {
  DataSet::DataSet() :
  Data::Basic::Set<Fasta::Entry>() {
  }

  DataSet::DataSet(const Specification::DataSet &s, bool revcomp, size_t n_seq) :
    Data::Basic::Set<Fasta::Entry>(s, revcomp, n_seq) { };

  DataSeries::DataSeries() :
  Data::Basic::Series<DataSet>() {
  }

  DataSeries::DataSeries(const std::string &name, const Specification::DataSets &paths, bool revcomp, size_t n_seq) :
  Data::Basic::Series<DataSet>(name, paths, revcomp, n_seq) {
  }

  DataCollection::DataCollection(const Specification::DataSets &paths, bool revcomp, size_t n_seq) :
  Data::Basic::Collection<DataSeries>(paths, revcomp, n_seq) {
  }
}

