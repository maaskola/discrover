#include <cmath>
#include <cstddef>
#include <map>

// this is based on python code from DREME

//  global constants
const double log10_ = log(10);
const double log_zero_ = -1e10;               // Zero on the log scale.
const double log_small_ = -0.5e10;            // Threshold below which everything is zero.
const double mm_nats_ = log(-log_zero_);      // dynamic range of double in NATs

//  Fisher's Exact Test
const double log0_99999999 = log(0.9999999);
const double log1_00000001 = log(1.00000001);

// Routines for computing the logarithm of a sum in log space.
double my_exp(double x) {
  if(x < log_small_)
    return 0.0;
  else
    return exp(x);
}

double log_sum1(double logx, double logy) {
  if((logx - logy) > mm_nats_)
    return logx;
  else
    return logx + log(1 + my_exp(logy - logx));
}

/*  Return the log(x+y) given log(x) and log(y). */
double log_sum(double logx, double logy) {
  if(logx > logy)
    return log_sum1(logx, logy);
  else
    return log_sum1(logy, logx);
}


// log gamma function using continued fractions
double lngamm(double z) {
  double x = 0.0;
  x = x + 0.1659470187408462e-06/(z+7.0);
  x = x + 0.9934937113930748e-05/(z+6.0);
  x = x - 0.1385710331296526    /(z+5.0);
  x = x + 12.50734324009056     /(z+4.0);
  x = x - 176.6150291498386     /(z+3.0);
  x = x + 771.3234287757674     /(z+2.0);
  x = x - 1259.139216722289     /(z+1.0);
  x = x + 676.5203681218835     /(z);
  x = x + 0.9999999999995183;
  return log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5);
}

// log n! computed using gamma function
std::map<size_t,double> lnfact_hash;
double lnfact(size_t n) {
  if(n<=1)
    return 0.0;

  auto iter = lnfact_hash.find(n);
  if(iter != end(lnfact_hash))
    return(iter->second);

  double result = lngamm(n+1.0);
  lnfact_hash[n] = result;
  return result;
}

// log binomial coefficient n choose k
double lnbico(size_t n, size_t k) {
  return lnfact(n)-lnfact(k)-lnfact(n-k);
}

double log_hyper_323(double n11, double n1_, double n_1, double n) {
  return lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1);
}


double _log_sprob = 0;
size_t _sn11, _sn1_, _sn_1, _sn;
double log_hyper0(size_t n11i, size_t n1_i, size_t n_1i, size_t ni) {
  // TODO
  // global _sn11, _sn1_, _sn_1, _sn, _log_sprob
  if(not ((n1_i|n_1i|ni)!=0)){
    if(not (n11i % 10 == 0)) {
      if(n11i==_sn11+1) {
        _log_sprob = _log_sprob + log( ((_sn1_-_sn11)/float(n11i))*((_sn_1-_sn11)/float(n11i+_sn-_sn1_-_sn_1)) );
        _sn11 = n11i;
        return _log_sprob;
      }
      if(n11i==_sn11-1) {
        _log_sprob = _log_sprob + log( ((_sn11)/float(_sn1_-n11i))*((_sn11+_sn-_sn1_-_sn_1)/float(_sn_1-n11i)) );
        _sn11 = n11i;
        return _log_sprob;
      }
    }
    _sn11 = n11i;
  } else {
    _sn11 = n11i;
    _sn1_=n1_i;
    _sn_1=n_1i;
    _sn=ni;
  }
  _log_sprob = log_hyper_323(_sn11,_sn1_,_sn_1,_sn);
  return _log_sprob;

  return 0;
}

double log_hyper(double n11) {
  return log_hyper0(n11,0,0,0);
}

 

struct Res {
  double prob;
  double sless;
  double sright;
  double sleft;
  double slarge;
};

/**  Computes Fisher's exact test based on a
  null-hypothesis distribution specified by the totals, and
  an observed distribution specified by b1 and b2, i.e.
  determines the probability of b's outcomes 1 and 2.

  Returns an immutable list consisting of the exact
  probability, and assorted p-values (sless, sright, sleft,
  slarg) based on the density.
*/
Res log_getFETprob(double a1, double a2, double b1, double b2) {
  double log_sless = log_zero_;
  double log_sright = log_zero_;
  double log_sleft = log_zero_;
  double log_slarge = log_zero_;
  double n = a1 + a2 + b1 + b2;

  double row1 = a1 + a2; // the row containing the null hypothesis
  double col1 = a1 + b1; // the column containing samples for outcome 1
  double mx = row1; // the maximum
  if(col1 < mx)
    mx = col1;
  double mi = row1 + col1 - n; // the minimum
  if(mi < 0)
    mi = 0;
  if(mi == mx)
    return {0, 0, 0, 0, 0};

  double log_prob = log_hyper0(a1, row1, col1, n);
  log_sleft = log_zero_;
  double log_p = log_hyper(mi);

  double i = mi + 1;
  while(log_p < (log0_99999999 + log_prob)) {
    log_sleft = log_sum(log_sleft, log_p);
    log_p = log_hyper(i);
    i = i + 1;
  }

  i = i - 1;
  if(log_p < (log1_00000001 + log_prob))
    log_sleft = log_sum(log_sleft, log_p);
  else
    i = i - 1;

  log_sright = log_zero_;
  log_p = log_hyper(mx);

  double j = mx - 1;
  while(log_p < (log0_99999999 + log_prob)) {
    log_sright = log_sum(log_sright, log_p);
    log_p = log_hyper(j);
    j = j - 1;
  }

  j = j + 1;
  if(log_p < (log1_00000001 + log_prob))
    log_sright = log_sum(log_sright, log_p);
  else
    j = j + 1;

  if(fabs(i - a1) < fabs(j - a1)) {
    log_sless = log_sleft;
    // log_slarge = (log_slarge, log_prob);
    log_slarge = log(1.0-exp(log_sleft));
    log_slarge = log_sum(log_slarge, log_prob);
  } else {
    // log_sless = log_sum(1.0, -log_sright);
    log_sless = log(1.0-exp(log_sright));
    // log_sless = (log_sless, log_prob);
    log_sless = log_sum(log_sless, log_prob);
    log_slarge = log_sright;
  }
  return {log_prob, log_sless, log_sright, log_sleft, log_slarge};
}

/*  Return log of hypergeometric pvalue of #pos >= p
    p = positive successes
    P = positives
    n = negative successes
    N = negatives
    log_pthresh = short-circuit if p will be greater
    */
double getLogFETPvalue(double p, double P, double n, double N, double log_pthresh) {
  // check that p-value is less than 0.5
  // if p/float(P) > n/float(N):
  // if (p * N > n * P):
  if((p * N > n * P) and log_hyper_323(p,P,n+p,N+P) < log_pthresh)
    // apply Fisher Exact test (hypergeometric p-value)
    return log_getFETprob(N-n, n, P-p, p).slarge;
  else
    return 0;          // pvalue = 1
}


