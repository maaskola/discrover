
#include "plasma_stats.hpp"
#include "../stats_config.hpp"
#include <cmath>

namespace Seeding {

}

using namespace std;

FisherExactTestResults fisher_exact_test(const matrix_t &X_,
    double hypothesizedOddsRatio,
    Alternative alternative,
    bool confInt,
    double confLevel,
    Verbosity verbosity) {
  FisherExactTestResults res;
  res.p_value = 1;
  res.log_p_value = 0;

  size_t nr = X_.size1();
  size_t nc = X_.size2();

  // check arguments
  if(nr < 2 or nc < 2) {
    cerr << "Fisher exact test: table must be at least 2x2." << endl;
    exit(-1);
  }
  count_matrix_t X(nr,nc);
  for(size_t i = 0; i < nr; i++)
    for(size_t j = 0; j < nr; j++) {
      double x = X_(i,j);
      if(x < 0 or std::isnan(x)) {
        cerr << "Fisher exact test: all entries of the table must be nonnegative and finite." << endl;
        exit(-1);
      }
      size_t y = x;
      double z = y;
      if(x != z) {
        cerr << "Fisher exact test: warning entry " << i << "," << j << " = " << x << " rounded to " << z << "." << endl;
      }
      X(i,j) = y;
    }

  if((nr == 2) and (nc == 2)) {
    if(not((finite(confLevel) != 0) and (confLevel > 0) and (confLevel < 1))) {
      cerr << "Fisher exact test: confidence level must be a number between 0 and 1." << endl;
      exit(-1);
    }
    if(std::isnan(hypothesizedOddsRatio) or (hypothesizedOddsRatio < 0)) {
      cerr << "Fisher exact test: hypothesized odds ratio must be a number not less than 0." << endl;
      exit(-1);
    }
  }

  if((nr != 2) or (nr != 2)) {
    cerr << "Fisher exact test: test not implemented for number or rows or columns > 2." << endl;
    exit(-1);
  }

  if((nr == 2) and (nr == 2)) {
    // if(hybrid)
    //   cerr << "Fisher exact test: warning - hybrid is ignored for 2x2 tables." << endl;

    double m = X(0,0) + X(1,0);
    double n = X(0,1) + X(1,1);
    double k = X(0,0) + X(0,1);
    double x = X(0,0);
    double lo = max<double>(0, k-n);
    double hi = min<double>(k, m);

    if(verbosity >= Verbosity::debug) {
      cerr << "Fisher exact test: m = " << m << endl;
      cerr << "Fisher exact test: n = " << n << endl;
      cerr << "Fisher exact test: k = " << k << endl;
      cerr << "Fisher exact test: x = " << x << endl;
      cerr << "Fisher exact test: lo = " << lo << endl;
      cerr << "Fisher exact test: hi = " << hi << endl;
    }

    const size_t l = hi - lo + 1;

    if(verbosity >= Verbosity::debug)
      cerr << "l = " << l << endl;
    vector<double> support(l);
    for(size_t i = 0; i < l; i++)
      support[i] = lo + i;

    vector<double> logdc(l);
    for(size_t i = 0; i < l; i++) {
      double z = dhyper(support[i], m, n, k, true);
      if(verbosity >= Verbosity::debug)
        cerr << "Fisher exact test: logcd = " << z << endl;
      logdc[i] = z;
    }

    auto dnhyper = [&](double ncp) {
      vector<double> d(l);
      const double lncp = log(ncp);
      for(size_t i = 0; i < l; i++)
        d[i] = logdc[i] + lncp * support[i];
      double max_d = *max_element(begin(d), end(d));
      double s = 0;
      for(auto &v: d)
        s += v = exp(v - max_d);
      for(auto &v: d)
        v /= s;
      return(d);
    };

    /*
    auto mnhyper = [&](double ncp) {
      vector<double> v;
      if(ncp == 0) {
        v.push_back(lo);
        return(v);
      }
      if(std::isinf(ncp) != 0) {
        v.push_back(hi);
        return(v);
      }
      v = dnhyper(ncp);
      for(size_t i = 0; i < v.size(); i++)
        v[i] = support[i] * v[i];
      return(v);
    };  */

    // phyper = function (q, m, n, k, lower.tail = TRUE, log.p = FALSE) 
 
    auto pnhyper = [&](double q, double ncp=1, bool upperTail=false) {
      if(verbosity >= Verbosity::debug)
        cerr << "pnhyper" << endl
          << q << " " << ncp << " " << upperTail << endl;
      if(ncp == 1) {
        if(upperTail)
          return phyper(x - 1, m, n, k, false, false);
        else
          return phyper(x, m, n, k, true, false);
      }
      if(ncp == 0) {
        if(upperTail)
          return double(q <= lo);
        else
          return double(q >= lo);
      }
      if(isinf(ncp == 1)) {
        if(upperTail)
          return double(q <= hi);
        else
          return double(q >= hi);
      }
      double s = 0;
      auto v = dnhyper(ncp);
      for(size_t i = 0; i < l; i++)
        if(upperTail) {
          if(support[i] >= q)
            s += v[i];
        } else {
          if(support[i] <= q)
            s += v[i];
        }
      return(s);
    };

    switch(alternative) {
      case Alternative::Less:
        res.p_value = pnhyper(x, hypothesizedOddsRatio);
        break;
      case Alternative::Greater:
        res.p_value = pnhyper(x, hypothesizedOddsRatio, true);
        break;
      case Alternative::TwoSided:
        if(hypothesizedOddsRatio == 0) {
          if(verbosity >= Verbosity::debug)
            cerr << "case x" << endl;
          res.p_value = x == lo;
        } else {
          if(std::isinf(hypothesizedOddsRatio) == 1) {
            if(verbosity >= Verbosity::debug)
              cerr << "case y" << endl;
            res.p_value = x == hi;
          } else {
            if(verbosity >= Verbosity::debug)
              cerr << "case z" << endl;
            double relErr = 1 + 1e-7;
            auto d = dnhyper(hypothesizedOddsRatio);
            if(verbosity >= Verbosity::debug) {
              cerr << "d =";
              for(auto &v: d)
                cerr << " " << v;
              cerr << endl;
            }
            res.p_value = 0;
            double e = d[x - lo] * relErr;
            if(verbosity >= Verbosity::debug)
              cerr << "e = " << e << endl;
            for(size_t i = 0; i < l; i++) {
              if(d[i] <= e)
                res.p_value += d[i];
            }
          }
        }
        break;
    }
  } else {
    cerr << "Fisher exact test: bla!" << endl;
  }

  res.log_p_value = log(res.p_value);

  if(verbosity >= Verbosity::debug) {
    cerr << "Fisher exact test: p-value = " << res.p_value << endl;
    cerr << "Fisher exact test: log p-value = " << res.log_p_value << endl;
  }

  return res;
}
