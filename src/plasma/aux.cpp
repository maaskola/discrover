
#include <iostream>
#include <cstring>
#include <boost/algorithm/string.hpp>
#include "aux.hpp"
#include "../sha1.h"

using namespace std;

string string_tolower (const string & str ) {
  string temp = str;
  range_tolower(begin(temp), end(temp));
  return(temp);
}

vector<string> tokenize(const string &s, const string &delim)
{
  vector<string> strs;
  boost::split(strs, s, boost::is_any_of(delim));
  return(strs);
}

vector<size_t> parse_list(const string &s)
{
  if(s.find_first_not_of("0123456789,-") != string::npos) {
    cout << "Please note that the format for the list specification only allows digits, '-', and ','." << endl;
    exit(-1);
  }
  vector<string> strs = tokenize(s, ",");

  vector<size_t> v;

  for(auto t : strs) {
    if(t.size() > 0) {
      vector<string> s2;
      boost::split(s2, t, boost::is_any_of("-"));
      size_t i = 0, j = 0;
      switch(s2.size()) {
        case 1:
          i = atoi(s2[0].c_str());
          v.push_back(i);
          break;
        case 2:
          i = atoi(s2[0].c_str());
          j = atoi(s2[1].c_str());
          for(size_t k = i; k <= j; k++)
            v.push_back(k);
          break;
        default:
          cout << "List format error: only one '-' is allowed in any group." << endl
            << "The offending group is '" << t << "'." << endl;
          exit(-1);
          break;
      }
    }
  }
  return(v);
}

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

