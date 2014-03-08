
#include <cmath>
#include <cassert>
#include <boost/math/distributions/chi_squared.hpp>
#include "association.hpp"
#include "../stats_config.hpp"

using namespace std;

const Verbosity association_verbosity = Verbosity::info;

inline vector_t row_sums(const matrix_t &m)
{
  size_t n = m.size1();
  vector_t v = zero_vector(n);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      v(i) += m(i,j);
  return(v);
}

inline vector_t col_sums(const matrix_t &m)
{
  size_t n = m.size2();
  vector_t v = zero_vector(n);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      v(j) += m(i,j);
  return(v);
}

double calc_mutual_information(const matrix_t &matrix, double pseudo_count, bool normalize, bool correction, bool variance)
{
  if(variance)
    correction = true;
  if(association_verbosity >= Verbosity::debug) {
    cout << "calc_mutual_information(matrix_t)" << endl;
    cout << "Calculating mutual information of this table:" << endl << matrix << endl;
  }
  matrix_t m = matrix;
  if(pseudo_count != 0) {
    m = m + pseudo_count;
    normalize = true;
  }
  double z = 1;
  if(normalize) {
    z = 0;
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        z += m(i,j);
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        m(i,j) /= z;
  }

  double mi = 0;
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      if(m(i,j) != 0)
        mi += m(i,j) * (log(m(i,j)) - log(rs(i)) - log(cs(j)));
  mi /= log(2.0);
  if(association_verbosity >= Verbosity::debug)
    cout << "mi = " << mi << endl;
  mi = max<double>(mi, 0);
  if(not correction)
    return(mi);
  double j = mi;
  if(not variance) {
    double corr = (m.size1() - 1) * (m.size2() - 1) / (z + 1) / 2;
    mi = mi + corr;
    return(mi);
  } else {
    double mi2 = 0;
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        if(m(i,j) != 0) {
          double x  = log(m(i,j)) - log(rs(i)) - log(cs(j));
          mi2 += m(i,j) * x * x;
        }
    mi2 /= log(2.0);

    double var = 1.0 / (z + 1) * (mi2 - j * j);
    return(var);
  }
}

double calc_g_test_from_mi(double mi, double n) {
  double g = 2 * n * log(2) * mi;
  return(g);
}

double calc_g_test(const matrix_t &matrix, double pseudo_count)
{
  if(association_verbosity >= Verbosity::debug) {
    cout << "calc_g_test(matrix_t)" << endl;
    cout << "Calculating g-test of this table:" << endl << matrix << endl;
  }
  size_t n = 0;
  for(size_t i = 0; i < matrix.size1(); i++)
    for(size_t j = 0; j < matrix.size2(); j++)
      n += matrix(i,j) + pseudo_count;
  double mi = calc_mutual_information(matrix, pseudo_count, true, false, false);
  double g = calc_g_test_from_mi(mi, n);
  return(g);
}

double corrected_pvalue(double score, double n, double df, double motif_len, Verbosity verbosity)
{
  double g = calc_g_test_from_mi(score, n);
  if(verbosity >= Verbosity::verbose)
    cout << "g = " << g << endl;
  double log_p = pchisq(g, df, false, true);
  if(verbosity >= Verbosity::verbose)
    cout << "log p(g) = " << log_p << endl;
  double cor_log_p = log(149) * motif_len + log_p;
  if(verbosity >= Verbosity::verbose)
    cout << "corrected log p(g) = " << cor_log_p << endl;
  return(cor_log_p);
}

double calc_log_likelihood_ratio(const matrix_t &matrix, double pseudo_count)
{
  if(association_verbosity >= Verbosity::debug) {
    cout << "calc_log_likelihood_ratio (matrix_t)" << endl;
    cout << "Calculating log likelihood ratio of this table:" << endl << matrix << endl;
  }
  size_t n = 0;
  for(size_t i = 0; i < matrix.size1(); i++)
    for(size_t j = 0; j < matrix.size2(); j++)
      n += matrix(i,j) + pseudo_count;
  double mi = calc_mutual_information(matrix, pseudo_count, true, false, false);
  double llr = - n * log(2) * mi;
  return(llr);
}

double chi_sq(const matrix_t &matrix, double pseudo_count=1.0) {
  matrix_t m = matrix;
  if(pseudo_count != 0) {
    for(size_t i = 0; i < m.size1(); i++)
      for(size_t j = 0; j < m.size2(); j++)
        m(i,j) += pseudo_count;
  }
  double z = 0;
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      z += m(i,j);

  double chisq = 0;
  vector_t rs = row_sums(m);
  vector_t cs = col_sums(m);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++) {
      double e = rs(i) * cs(j) / z;
      // chisq += (m(i,j) - e) * (m(i,j) - e) / m(i,j);
      chisq += (m(i,j) - e) * (m(i,j) - e) / e;
    }
  return(chisq);
}

/** Computes mutual information of a contingency table
 * a,b,c,d are assumed to be non-negative
 */
double calc_mutual_information(double A, double B, double C, double D, bool normalize)
{
  if(association_verbosity >= Verbosity::debug)
    cout << "calc_mutual_information(A, B, C, D)" << endl;
  double a = A, b = B, c = C, d = D;
  if(normalize) {
    double n = a + b + c + d;
    a /= n;
    b /= n;
    c /= n;
    d /= n;
  }
  double mi = 0;
  double lab = log(a+b);
  double lac = log(a+c);
  double lbd = log(b+d);
  double lcd = log(c+d);
  mi += a * (log(a) - lab - lac);
  mi += b * (log(b) - lab - lbd);
  mi += c * (log(c) - lac - lcd);
  mi += d * (log(d) - lbd - lcd);
  mi /= log(2);
  if(association_verbosity >= Verbosity::debug)
    cout << "mi = " << mi << endl;
  return(mi);
}

/** Computes mutual information of a contingency table
 * a,b,c,d are assumed to be non-negative
 */
double calc_matthews_correlation_coefficient(double tp, double fp, double fn, double tn)
{
  // cout << "mcc(" << tp << ", " << fp << ", " << fn << ", " << tn << ")" << endl;
  double a = tp+fp;
  double b = tp+fn;
  double c = tn+fp;
  double d = tn+fn;
  if(a == 0 or b == 0 or c == 0 or d == 0)
    return 0;
  double mcc = (tp * tn - fp * fn) / sqrt(a*b*c*d);
  return(mcc);
}

double calc_matthews_correlation_coefficient(const confusion_matrix &m)
{
  return(calc_matthews_correlation_coefficient(m.true_positives, m.false_negatives, m.false_positives, m.true_negatives));
}

double calc_mutual_information(const confusion_matrix &m, bool normalize)
{
  if(association_verbosity >= Verbosity::debug)
    cout << "calc_mutual_information(confusion_matrix)" << endl;
  double mi = calc_mutual_information(m.true_positives, m.false_negatives, m.false_positives, m.true_negatives, normalize);
  if(association_verbosity >= Verbosity::debug)
    cout << "mi = " << mi << endl;
  return(mi);
}

double calc_rank_information(vector_t posterior, double pseudo_count)
{
  vector_t p = posterior;
  vector_t q = posterior;
  size_t n = posterior.size();

  // forward cumulative sum
  double z = 0;
  for(size_t i = 0; i < n; i++)
    z = p[i] += z;

  // backward cumulative sum
  z = 0;
  for(size_t i = 0; i < n; i++)
    z = q[n - 1 - i] += z;

  // shift back one position and set the last to zero
  for(size_t i = 0; i < n-1; i++)
    q[i] = q[i+1];
  q[n-1] = 0;

  vector_t r = p;
  vector_t s = q;
  for(size_t i = 0; i < n; i++)
    r[i] = i + 1 - r[i];
  for(size_t i = 0; i < n; i++)
    s[i] = n - 1 - i - s[i];

  double ric = 0;
  for(size_t i = 0; i < n-1; i++) {
    double mi = calc_mutual_information(p[i] + pseudo_count, q[i] + pseudo_count, r[i] + pseudo_count, s[i] + pseudo_count, true);
    ric += mi;
    // cout << "debug " << i << " " << p[i] << " " << q[i] << " " << r[i] << " " << s[i] << " " << p[i] + q[i] + r[i] + s[i] << " " << mi << " " << ric << endl;
  }
  // cout << "debug ric = " <<  ric << endl;
  ric /= n;
  return(ric);
}



/*
double OccurrenceStatistics::chisq_p(double x, size_t df) const
{
  using boost::math::chi_squared;
  using boost::math::quantile;
  using boost::math::complement;
  chi_squared dist(df);
  return(cdf(complement(dist, x)));
}
*/

