/* =====================================================================================
 * Copyright (c) 2011, Jonas Maaskola
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 *
 *       Filename:  hmm_linesearch.cpp
 *
 *    Description:  Line searching routines for gradient learning with HMMs
 *
 *        Created:  Wed Oct 8 21:06:55 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <iomanip>
#include "../timer.hpp"
#include "../aux.hpp"
#include "logistic.hpp"
#include "../matrix_inverse.hpp"
#include "polyfit.hpp"
#include "hmm.hpp"

using namespace std;

const double p5 = 0.5;
const double p66 = 0.66;

template <class X>
X max3(const X &a, const X &b, const X &c) {
  return max<X>(a, max<X>(b, c));
};

double update(double current, double step) {
  double x = logit(current);  // log(current) - log(1 - current);
  // cout << "setting: current " << current << " x = " << x << " e(x) = " <<
  // (exp(x) / (exp(x)+1)) << endl;
  double y = x + step;
  double z = exp(y);
  double u = z / (1 + z);
  // cout << "setting: u = " << u << endl;
  // return current;
  return u;
}

bool mcstep(double &stx, double &fx, double &dx, double &sty, double &fy,
            double &dy, double &stp, double fp, double dp, bool &brackt,
            double stpmin, double stpmax, Verbosity verbosity);

HMM HMM::build_trial_model(const Gradient &gradient, double alpha,
                           const Training::Task &task) const {
  matrix_t t = transition;
  matrix_t e = emission;

  log_transform(t);
  log_transform(e);

  double z = sqrt(scalar_product(gradient, gradient));

  matrix_t t_step = t;
  if (not task.targets.transition.empty())
    t_step += alpha * gradient.transition / z;

  matrix_t e_step = e;
  if (not task.targets.emission.empty())
    e_step += alpha * gradient.emission / z;

  // TODO idea: center on 0 before exp_transforming
  // note: in computing the mean (that will be subtracted for centering) skip
  // those that are negative infinite
  exp_transform(t_step);
  exp_transform(e_step);

  // fix excessive values and normalize
  for (size_t i = 0; i < t.size1(); i++) {
    // fix excessive transition probability values to something reasonable
    for (size_t j = 0; j < t.size2(); j++)
      if (t_step(i, j) > numeric_limits<double>::max() / 100)
        t_step(i, j) = numeric_limits<double>::max() / 100;

    // normalize transition probabilities
    normalize_transition(t_step);

    // fix excessive emission probability values to something reasonable
    for (size_t j = 0; j < e.size2(); j++)
      if (e_step(i, j) > numeric_limits<double>::max() / 100)
        e_step(i, j) = numeric_limits<double>::max() / 100;

    // normalize emission probabilities
    normalize_emission(e_step);
  }

  if (verbosity >= Verbosity::verbose) {
    cerr << "Normalized transition step: " << t_step << endl;
    cerr << "Normalized emission step: " << e_step << endl;
  }

  HMM trial_hmm(*this);
  if (not task.targets.transition.empty())
    trial_hmm.transition = t_step;
  if (not task.targets.emission.empty())
    trial_hmm.emission = e_step;

  return trial_hmm;
}

inline double psi(double v1, double v2, const Gradient &g, double alpha,
                  double mu, const Gradient &direction, Verbosity verbosity) {
  double d = dderiv(direction, g);
  double psi = v2 - v1 - alpha * mu * d;
  if (verbosity >= Verbosity::debug)
    cout << "alpha = " << alpha << " v1 = " << v1 << " v2 = " << v2
         << " mu = " << mu << " d = " << d << " psi = " << psi << endl;
  return psi;
}

inline double psi_gradient(const Gradient &g1, const Gradient &g2, double alpha,
                           double mu, const Gradient &direction,
                           Verbosity verbosity) {
  double d1 = dderiv(direction, g1);
  double d2 = dderiv(direction, g2);
  double psi = d2 - mu * d1;
  if (verbosity >= Verbosity::debug)
    cout << "alpha = " << alpha << " mu = " << mu << " d1 = " << d1
         << " d2 = " << d2 << " psi gradient = " << psi << endl;
  return psi;
}

bool sufficient_increase_condition(double v1, double v2, const Gradient &g,
                                   double alpha, double mu,
                                   const Gradient &direction,
                                   Verbosity verbosity) {
  bool criterion = psi(v1, v2, g, alpha, mu, direction, verbosity) >= 0;
  if (verbosity >= Verbosity::debug)
    cout << "Sufficient increase criterion is "
         << (criterion ? "FULFILLED" : "NOT FULFILLED") << endl;
  return criterion;
}

bool curvature_condition(double v1, double v2, const Gradient &g1,
                         const Gradient &g2, double alpha, double eta,
                         const Gradient &direction, Verbosity verbosity) {
  double d1 = dderiv(direction, g1);
  double d2 = dderiv(direction, g2);
  bool criterion = fabs(d2) <= eta * fabs(d1);
  if (verbosity >= Verbosity::debug)
    cout << "Curvature criterion is "
         << (criterion ? "FULFILLED" : "NOT FULFILLED") << ": "
         << "|" << d2 << "| " << (criterion ? "<=" : ">") << " " << eta
         << " * |" << d1 << "| = " << (eta * fabs(d1)) << endl;
  return criterion;
}

Gradient normalize(const Gradient &gradient) {
  Gradient u = gradient;
  double norm = sqrt(scalar_product(u, u));
  u.transition /= norm;
  u.emission /= norm;
  return u;
}

/** Line-searching algorithm due to
 * Jorge J. Moré and David J. Thuente.
 * Line Search Algorithms with Guaranteed Sufficient Decrease
 * ACM Transactions on Mathematical Software (TOMS)
 * Volume 20 Issue 3, Sept. 1994
 * Pages 286-307
 * 
 * INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
 *
 *        INFO = 0  IMPROPER INPUT PARAMETERS.
 *
 *        INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
 *
 *        INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
 *                  DIRECTIONAL DERIVATIVE CONDITION HOLD.
 *
 *        INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
 *                  IS AT MOST XTOL.
 *
 *        INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
 *
 *        INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
 *
 *        INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
 *
 *        INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
 *                  THERE MAY NOT BE A STEP WHICH SATISFIES THE
 *                  SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
 *                  TOLERANCES MAY BE TOO SMALL.
 *
 **/
pair<double, HMM> HMM::line_search_more_thuente(
    const Data::Collection &collection, const Gradient &initial_gradient_,
    double initial_score, int &info, const Training::Task &task,
    const Options::HMM &options) const {
  const Verbosity verbo = verbosity;
  // const Verbosity verbo = Verbosity::verbose;
  // const Verbosity verbo = Verbosity::debug;

  info = 0;
  bool failed = false;

  Gradient initial_gradient = initial_gradient_;

  // check option validity
  if (options.line_search.max_steps == 0) {
    cout << "Error: line search option max_steps must be a positive integer."
         << endl;
    failed = true;
  }
  if (options.line_search.mu < 0) {
    cout << "Error: line search option mu must be non-negative." << endl;
    failed = true;
  }
  if (options.line_search.eta < 0) {
    cout << "Error: line search option eta must be non-negative." << endl;
    failed = true;
  }
  if (options.line_search.eta <= options.line_search.mu) {
    cout << "Error: line search option eta must be greater than mu." << endl;
    failed = true;
  }

  size_t nfev = 1;

  const Gradient direction = normalize(initial_gradient);

  double dginit = dderiv(direction, initial_gradient);
  if (dginit <= 0) {
    cout << "Direction is not an ascent direction!" << endl;
    failed = true;
  }
  if (failed)
    return pair<double, HMM>(initial_score, *this);

  // TODO determine alpha_min and alpha_max
  const double min_score = initial_score
                           / (1 - options.termination.delta_tolerance);
  const double stpmin = (min_score - initial_score)
                        / dderiv(direction, initial_gradient);
  // const double stpmin = (min_score - initial_score) / options.line_search.mu / dderiv(direction, initial_gradient);
  if (verbo >= Verbosity::debug) {
    cout << "min_score = " << min_score << endl;
    cout << "stpmin= " << stpmin << endl;
  }
  // double alpha_max = 1.0 / options.line_search.mu * (min_score - initial_score) / gradient_sum(initial_gradient);

  double upper_limit = 10000;
  if (task.measure == Measure::MutualInformation
      or task.measure == Measure::MatthewsCorrelationCoefficient)
    upper_limit = (1 - initial_score) / options.line_search.mu
                  / dderiv(direction,
                           initial_gradient);  // Assume the maximal score is 1
  if (verbo >= Verbosity::debug)
    cout << "Upper limit = " << upper_limit << endl;
  double strict_upper_limit = min(upper_limit + 1, upper_limit * 1.1);
  if (verbo >= Verbosity::debug)
    cout << "Strict upper limit = " << strict_upper_limit << endl;

  double stpmax = strict_upper_limit;

  const double xtrapf = 4;
  const double xtol = 1e-10;
  const double ftol = options.line_search.mu;
  const double gtol = options.line_search.eta;
  const size_t maxfev = options.line_search.max_steps;

  double finit = initial_score;
  double dgtest = ftol * dginit;
  double width = stpmax - stpmin;
  double width1 = width / p5;

  if (verbo >= Verbosity::debug) {
    cout << "ftol = " << ftol << " gtol = " << gtol << endl;
    cout << "dginit = " << dginit << " dgtest = " << dgtest << endl;
  }

  // THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
  // FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
  // THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
  // FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
  // THE INTERVAL OF UNCERTAINTY.
  // THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
  // FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.

  double stp = 1;
  double stx = 0;
  double fx = initial_score;
  double dgx = dginit;
  double sty = 0;
  double fy = initial_score;
  double dgy = dginit;

  bool brackt = false;
  bool stage1 = true;

  // START OF ITERATIONS
  while (true) {
    if (verbo >= Verbosity::verbose) {
      cout << "Starting line search iteration. stp = " << stp << endl
           << "stx = " << stx << " fx = " << fx << " dgx = " << dgx << endl;
      if (brackt)
        cout << "sty = " << sty << " fy = " << fy << " dgy = " << dgy << endl;
    }

    // SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
    // TO THE PRESENT INTERVAL OF UNCERTAINTY.
    double stmin, stmax;
    if (brackt) {
      stmin = min(stx, sty);
      stmax = max(stx, sty);
    } else {
      stmin = stx;
      stmax = stp + xtrapf * (stp - stx);
    }

    if (verbo >= Verbosity::verbose)
      cout << "stmin = " << stmin << " stmax = " << stmax << endl;

    // FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
    stp = max(stp, stpmin);
    stp = min(stp, stpmax);

    // IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
    // STP BE THE LOWEST POINT OBTAINED SO FAR.
    if ((brackt and (stp <= stmin or stp >= stmax)) or nfev >= maxfev - 1
        or (brackt and stmax - stmin <= xtol * stmax))
      stp = stx;

    // EVALUATE THE FUNCTION AND GRADIENT AT STP
    // AND COMPUTE THE DIRECTIONAL DERIVATIVE.
    // We return to main program to obtain F and G.
    nfev++;
    double f;
    HMM trial_hmm = build_trial_model(initial_gradient, stp, task);
    Gradient trial_gradient
        = trial_hmm.compute_gradient(collection, f, task, options.weighting);
    if (verbo >= Verbosity::verbose)
      cout << "Function and gradient evaluation!" << endl;

    double dg = dderiv(direction, trial_gradient);
    double ftest1 = finit + stp * dgtest;

    if (verbo >= Verbosity::verbose)
      cout << "stp = " << stp << " f = " << f << " dg = " << dg << endl;
    if (verbo >= Verbosity::debug)
      cout << "ftest1 = " << ftest1 << endl;

    // TEST FOR CONVERGENCE.
    if ((brackt and (stp <= stmin or stp >= stmax)))
      info = 6;
    if (stp == stpmax and f >= ftest1 and dg >= dgtest)
      info = 5;
    if (stp == stpmin and (f < ftest1 or dg <= dgtest))
      info = 4;
    if (nfev >= maxfev)
      info = 3;
    if (brackt and stmax - stmin <= xtol * stmax)
      info = 2;
    if (f >= ftest1 and fabs(dg) <= gtol * fabs(dginit))
      info = 1;

    // CHECK FOR TERMINATION.
    if (info != 0) {
      if (verbo >= Verbosity::info)
        cout << "Fnc & grad evaluations in line search          " << nfev
             << endl;
      return pair<double, HMM>(f, trial_hmm);
    }

    // IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
    // FUNCTION HAS A NONNEGATIVE VALUE AND NONPOSITIVE DERIVATIVE.
    if (stage1 and f >= ftest1 and dg <= min(ftol, gtol) * dginit)
      stage1 = false;

    // A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
    // WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
    // FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
    // DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
    // OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.

    if (stage1 and f >= fx and f < ftest1) {
      // DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
      if (verbo >= Verbosity::verbose)
        cout << "Using modified function and derivative values." << endl;
      double fm = f - stp * dgtest;
      double fxm = fx - stx * dgtest;
      double fym = fy - sty * dgtest;
      double dgm = dg - dgtest;
      double dgxm = dgx - dgtest;
      double dgym = dgy - dgtest;

      // CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
      // AND TO COMPUTE THE NEW STEP.
      if (not mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt,
                     stmin, stmax, verbo))
        // the intial parameters of mcstep were wrong
        return pair<double, HMM>(initial_score, *this);

      // RESET THE FUNCTION AND GRADIENT VALUES FOR F.
      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;

    } else {
      // CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
      // AND TO COMPUTE THE NEW STEP.
      if (verbo >= Verbosity::verbose)
        cout << "Using original function and derivative values." << endl;
      if (not mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin,
                     stmax, verbo))
        // the intial parameters of mcstep were wrong
        return pair<double, HMM>(initial_score, *this);
    }

    // FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
    // INTERVAL OF UNCERTAINTY.
    if (brackt) {
      if (fabs(sty - stx) >= p66 * width1)
        stp = stx + p5 * (sty - stx);
      width1 = width;
      width = fabs(sty - stx);
    }
  }
}

bool mcstep(double &stx, double &fx, double &dx, double &sty, double &fy,
            double &dy, double &stp, double fp, double dp, bool &brackt,
            double stpmin, double stpmax, Verbosity verbosity) {
  if (verbosity >= Verbosity::debug)
    cout << "MCSTEP " << stx << " " << fx << " " << dx << " " << sty << " "
         << fy << " " << dy << " " << stp << " " << fp << " " << dp << " "
         << brackt << " " << stpmin << " " << stpmax << endl;
  bool bound = false;
  // SUBROUTINE MCSTEP
  //
  // THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
  // A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
  // A MAXIMIZER OF THE FUNCTION.
  //
  // THE PARAMETER STX CONTAINS THE STEP WITH THE GREATEST FUNCTION
  // VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
  // ASSUMED THAT THE DERIVATIVE AT STX IS POSITIVE IN THE
  // DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
  // MAXIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
  // WITH ENDPOINTS STX AND STY.
  //
  // THE SUBROUTINE STATEMENT IS
  //
  //   SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
  //                    STPMIN,STPMAX,INFO)
  //
  // WHERE
  //
  //   STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
  //     THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
  //     SO FAR. THE DERIVATIVE MUST BE POSITIVE IN THE DIRECTION
  //     OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE EQUAL
  //     SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
  //
  //   STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
  //     THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
  //     THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
  //     UPDATED APPROPRIATELY.
  //
  //   STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
  //     THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
  //     IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
  //     BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
  //
  //   BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MAXIMIZER
  //     HAS BEEN BRACKETED. IF THE MAXIMIZER HAS NOT BEEN BRACKETED
  //     THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MAXIMIZER
  //     IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
  //
  //   STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
  //     AND UPPER BOUNDS FOR THE STEP.
  //
  //   INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
  //     IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
  //     ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
  //     INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
  //     
  // ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
  // JORGE J. MORE', DAVID J. THUENTE

  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  // CHECK THE INPUT PARAMETERS FOR ERRORS.
  if ((brackt and (stp <= min(stx, sty) or stp >= max(stx, sty)))
      or dx * (stp - stx) <= 0 or stpmax < stpmin) {
    cout << "Error: Input values are not valid." << endl;
    cout << "brackt = " << brackt << endl;
    cout << "stp = " << stp << endl;
    cout << "stx = " << stx << endl;
    cout << "sty = " << sty << endl;
    cout << "dx = " << dx << endl;
    cout << "stpmax = " << stpmax << endl;
    cout << "stpmin = " << stpmin << endl;
    return false;
  }

  // DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
  sgnd = dp * (dx / fabs(dx));

  // FIRST CASE. A LOWER FUNCTION VALUE.
  // THE MAXIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
  // TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
  // ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
  if (fp < fx) {
    if (verbosity >= Verbosity::verbose)
      cout << "FIRST CASE. A LOWER FUNCTION VALUE." << endl;
    bound = true;
    stpc = cubic_extrema(interpolate3(stx, stp, fx, fp, dx, dp),
                         Verbosity::verbose).second;
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
    if (verbosity >= Verbosity::debug)
      cout << "stpc = " << stpc << " stpq = " << stpq << endl;
    if (fabs(stpc - stx) < fabs(stpq - stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / 2;
    brackt = true;
  } else if (sgnd < 0) {
    // SECOND CASE. A HIGHER FUNCTION VALUE AND DERIVATIVES OF
    // OPPOSITE SIGN. THE MAXIMUM IS BRACKETED. IF THE CUBIC
    // STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
    // THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
    if (verbosity >= Verbosity::verbose)
      cout << "SECOND CASE. A HIGHER FUNCTION VALUE AND DERIVATIVES OF "
              "OPPOSITE SIGN." << endl;
    bound = false;
    stpc = cubic_extrema(interpolate3(stx, stp, fx, fp, dx, dp),
                         Verbosity::verbose).second;
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (verbosity >= Verbosity::debug)
      cout << "stpc = " << stpc << " stpq = " << stpq << endl;
    if (fabs(stpc - stp) > fabs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  } else if (fabs(dp) < fabs(dx)) {
    // THIRD CASE. A HIGHER FUNCTION VALUE, DERIVATIVES OF THE
    // SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
    // THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
    // IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
    // IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
    // EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
    // COMPUTED AND IF THE MAXIMUM IS BRACKETED THEN THE THE STEP
    // CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
    if (verbosity >= Verbosity::verbose)
      cout << "THIRD CASE. A HIGHER FUNCTION VALUE, DERIVATIVES OF THE SAME "
              "SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES." << endl;
    bound = true;
    // theta = 3*(fx - fp)/(stp - stx) + dx + dp;
    // s = max3(fabs(theta),fabs(dx),fabs(dp));
    // // THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
    // // TO INFINITY IN THE DIRECTION OF THE STEP.
    // gamma = s*sqrt(max(0.0,pow(theta/s,2) - (dx/s)*(dp/s)));
    // cout << "theta = " << theta << " s = " << s << " gamma = " << gamma <<
    // endl;
    // if(stp > stx) gamma = -gamma;
    // p = (gamma - dp) + theta;
    // q = (gamma + (dx - dp)) + gamma;
    // r = p/q;
    // cout << "gamma = " << gamma << " p = " << p << " q = " << q << " r = " <<
    // r << endl;
    double a_c = cubic_extrema(interpolate3(stx, stp, fx, fp, dx, dp),
                               Verbosity::verbose).second;
    // if(r < 0.0 and gamma != 0.0)
    //   stpc = stp + r*(stx - stp);
    if (fabs(a_c - stx) > fabs(stp - stx))
      stpc = a_c;
    else if (stp > stx)
      stpc = stpmax;
    else
      stpc = stpmin;
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (verbosity >= Verbosity::debug)
      cout << "stpc = " << stpc << " stpq = " << stpq << endl;
    if (brackt) {
      if (fabs(stp - stpc) < fabs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    } else {
      if (fabs(stp - stpc) > fabs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    }
  } else {
    // FOURTH CASE. A HIGHER FUNCTION VALUE, DERIVATIVES OF THE
    // SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
    // NOT DECREASE. IF THE MAXIMUM IS NOT BRACKETED, THE STEP
    // IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
    if (verbosity >= Verbosity::verbose)
      cout << "FOURTH CASE. A HIGHER FUNCTION VALUE, DERIVATIVES OF THE SAME "
              "SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES NOT DECREASE."
           << endl;
    bound = false;
    if (brackt) {
      theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      s = max3(fabs(theta), fabs(dy), fabs(dp));
      gamma = s * sqrt(pow(theta / s, 2) - (dy / s) * (dp / s));
      if (stp > sty)
        gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      if (verbosity >= Verbosity::debug)
        cout << "stpc = " << stpc << endl;
      stpf = stpc;
    } else if (stp > stx)
      stpf = stpmax;
    else
      stpf = stpmin;
  }
  if (verbosity >= Verbosity::debug)
    cout << "Computed next step." << endl << "stpf = " << stpf << endl;

  // UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
  // DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
  if (fp < fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  if (verbosity >= Verbosity::verbose)
    cout << "Updated interval of uncertainty." << endl << "stx = " << stx
         << " fx = " << fx << " dx = " << dx << endl << "sty = " << sty
         << " fy = " << fy << " dy = " << dy << endl;

  // COMPUTE THE NEW STEP AND SAFEGUARD IT.
  stpf = min(stpmax, stpf);
  stpf = max(stpmin, stpf);
  stp = stpf;
  if (brackt and bound) {
    if (sty > stx)
      stp = min(stx + p66 * (sty - stx), stp);
    else
      stp = max(stx + p66 * (sty - stx), stp);
  }

  if (verbosity >= Verbosity::verbose)
    cout << "Safeguarded step." << endl << "stp = " << stp << endl;
  return true;
}
