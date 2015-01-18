
#include "chisq.hpp"
#include "pgamma.hpp"
#include <iostream>

using namespace std;

int main(int argc, const char** argv) {
  cout << "df\tx\tdp\tlogdp\tp\tlogp" << endl;
  for (size_t df = 0; df < 20; df++)
    for (size_t i = 0; i < 20; i++) {
      double x = i;
      double dp = dchisq(x, df, false);
      double logdp = dchisq(x, df, true);
      double p = pchisq(x, df, false, false);
      double logp = pchisq(x, df, false, true);
      cout << df << "\t" << i << "\t" << dp << "\t" << logdp << "\t" << p
           << "\t" << logp << endl;
    }
}
