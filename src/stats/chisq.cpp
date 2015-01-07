
#include "chisq.hpp"
// #include "pgamma.cpp"
#include <limits>
#include <cmath>
#include <cfloat>
#include <boost/math/special_functions/gamma.hpp>

/*
 *  stirlerr
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */

double stirlerr(double n) {
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2*pi)) */

#define S0 0.083333333333333333333        /* 1/12 */
#define S1 0.00277777777777777777778      /* 1/360 */
#define S2 0.00079365079365079365079365   /* 1/1260 */
#define S3 0.000595238095238095238095238  /* 1/1680 */
#define S4 0.0008417508417508417508417508 /* 1/1188 */

  /*
     error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
     */
  const static double sferr_halves[31] = {
      0.0,                           /* n=0 - wrong, place holder only */
      0.1534264097200273452913848,   /* 0.5 */
      0.0810614667953272582196702,   /* 1.0 */
      0.0548141210519176538961390,   /* 1.5 */
      0.0413406959554092940938221,   /* 2.0 */
      0.03316287351993628748511048,  /* 2.5 */
      0.02767792568499833914878929,  /* 3.0 */
      0.02374616365629749597132920,  /* 3.5 */
      0.02079067210376509311152277,  /* 4.0 */
      0.01848845053267318523077934,  /* 4.5 */
      0.01664469118982119216319487,  /* 5.0 */
      0.01513497322191737887351255,  /* 5.5 */
      0.01387612882307074799874573,  /* 6.0 */
      0.01281046524292022692424986,  /* 6.5 */
      0.01189670994589177009505572,  /* 7.0 */
      0.01110455975820691732662991,  /* 7.5 */
      0.010411265261972096497478567, /* 8.0 */
      0.009799416126158803298389475, /* 8.5 */
      0.009255462182712732917728637, /* 9.0 */
      0.008768700134139385462952823, /* 9.5 */
      0.008330563433362871256469318, /* 10.0 */
      0.007934114564314020547248100, /* 10.5 */
      0.007573675487951840794972024, /* 11.0 */
      0.007244554301320383179543912, /* 11.5 */
      0.006942840107209529865664152, /* 12.0 */
      0.006665247032707682442354394, /* 12.5 */
      0.006408994188004207068439631, /* 13.0 */
      0.006171712263039457647532867, /* 13.5 */
      0.005951370112758847735624416, /* 14.0 */
      0.005746216513010115682023589, /* 14.5 */
      0.005554733551962801371038690  /* 15.0 */
  };
  double nn;

  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn)
      return sferr_halves[(int)nn];
    return boost::math::lgamma(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI;
  }

  nn = n * n;
  if (n > 500)
    return (S0 - S1 / nn) / n;
  if (n > 80)
    return (S0 - (S1 - S2 / nn) / nn) / n;
  if (n > 35)
    return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  /* 15 < n <= 35 : */
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

double bd0(double x, double np) {
  double ej, s, s1, v;
  int j;

  if (!std::isfinite(x) || !std::isfinite(np) || np == 0.0)
    return NAN;
  // if(!std::isfinite(x) || !std::isfinite(np) || np == 0.0) ML_ERR_return_NAN;

  if (fabs(x - np) < 0.1 * (x + np)) {
    v = (x - np) / (x + np);
    s = (x - np) * v; /* s using v -- change by MM */
    ej = 2 * x * v;
    v = v * v;
    for (j = 1;; j++) { /* Taylor series */
      ej *= v;
      s1 = s + ej / ((j << 1) + 1);
      if (s1 == s) /* last term was effectively 0 */
        return s1;
      s = s1;
    }
  }
  /* else:  | x - np |  is not too small */
  return x * log(x / np) + np - x;
}

double dpois_raw(double x, double lambda, int give_log) {
  /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
           lambda >= 0
           */
  if (lambda == 0)
    return (x == 0) ? (give_log ? 0 : 1)
                    : (give_log ? -std::numeric_limits<double>::infinity() : 0);
  if (!std::isfinite(lambda))
    return give_log ? -std::numeric_limits<double>::infinity() : 0;
  if (x < 0)
    return give_log ? -std::numeric_limits<double>::infinity() : 0;
  if (x <= lambda * DBL_MIN)
    return give_log ? -lambda : exp(-lambda);
  if (lambda < x * DBL_MIN)
    return give_log ? -lambda + x * log(lambda) - boost::math::lgamma(x + 1) :
                    //-lambda + x*log(lambda) -lgammafn(x+1) :
               exp(-lambda + x * log(lambda) - boost::math::lgamma(x + 1));
  // exp(-lambda + x*log(lambda) -lgammafn(x+1)));
  return give_log ? -0.5 * log(2 * M_PI * x) - stirlerr(x) - bd0(x, lambda)
                   : exp(-stirlerr(x) - bd0(x, lambda)) / sqrt(2 * M_PI * x);
}

double dgamma(double x, double shape, double scale, int give_log) {
  double pr;
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
    return x + shape + scale;
#endif
  if (shape < 0 || scale <= 0)
    return NAN;
  // if (shape < 0 || scale <= 0) ML_ERR_return_NAN;
  if (x < 0)
    return give_log ? -std::numeric_limits<double>::infinity() : 0;
  if (shape == 0) /* point mass at 0 */
    return (x == 0) ? std::numeric_limits<double>::infinity()
                    : (give_log ? -std::numeric_limits<double>::infinity() : 0);
  if (x == 0) {
    if (shape < 1)
      return std::numeric_limits<double>::infinity();
    if (shape > 1)
      return give_log ? -std::numeric_limits<double>::infinity() : 0;
    /* else */
    return give_log ? -log(scale) : 1 / scale;
  }

  if (shape < 1) {
    pr = dpois_raw(shape, x / scale, give_log);
    return give_log ? pr + log(shape / x) : pr * shape / x;
  }
  /* else  shape >= 1 */
  pr = dpois_raw(shape - 1, x / scale, give_log);
  return give_log ? pr - log(scale) : pr / scale;
}

double dchisq(double x, double df, int give_log) {
  return dgamma(x, df / 2., 2., give_log);
}
