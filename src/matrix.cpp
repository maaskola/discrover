#include <fstream>
#include "aux.hpp"

using namespace std;

vector_t rowsums(const matrix_t &m) {
  size_t n1 = m.size1();
  size_t n2 = m.size2();
  vector_t v(n1);
  for (size_t i = 0; i < n1; i++) {
    v(i) = 0;
    for (size_t j = 0; j < n2; j++)
      v(i) += m(i, j);
  }
  return v;
}

vector_t colsums(const matrix_t &m) {
  size_t n1 = m.size2();
  size_t n2 = m.size1();
  vector_t v(n1);
  for (size_t i = 0; i < n1; i++) {
    v(i) = 0;
    for (size_t j = 0; j < n2; j++)
      v(i) += m(j, i);
  }
  return v;
}

matrix_t qki_consensus_pssm() {
  matrix_t p = zero_matrix(8, 4);
  p(0, 3) = 1;
  p(1, 0) = 1;
  p(2, 2) = 1;
  p(3, 3) = 1;
  p(4, 0) = 1;
  p(5, 0) = 1;
  p(6, 2) = 1;
  p(7, 0) = 1;
  p = p + scalar_matrix(8, 4, 0.01);
  //  p = p / 100;
  return p;
}

matrix_t pum_consensus_pssm() {
  matrix_t p = zero_matrix(8, 4);
  p(0, 3) = 1;
  p(1, 2) = 1;
  p(2, 3) = 1;
  p(3, 0) = 1;
  p(4, 0) = 0.5;
  p(4, 3) = 0.5;
  p(5, 0) = 1;
  p(6, 3) = 1;
  p(7, 0) = 1;
  p = p + scalar_matrix(8, 4, 0.01);
  //  p = p / 100;
  return p;
}

matrix_t operator+(const matrix_t &a, double x) {
  matrix_t b(a);
  for (size_t i = 0; i < b.size1(); i++)
    for (size_t j = 0; j < b.size2(); j++)
      b(i, j) = a(i, j) + x;
  return b;
}

matrix_t operator-(const matrix_t &a, double x) {
  matrix_t b(a);
  for (size_t i = 0; i < b.size1(); i++)
    for (size_t j = 0; j < b.size2(); j++)
      b(i, j) = a(i, j) - x;
  return b;
}

void exp_transform(double &x) { x = exp(x); }

void exp_transform(matrix_t &m) {
  for (size_t i = 0; i < m.size1(); i++)
    for (size_t j = 0; j < m.size2(); j++)
      m(i, j) = exp(m(i, j));
}

void log_transform(double &x) { x = log(x); }

void log_transform(matrix_t &m) {
  for (size_t i = 0; i < m.size1(); i++)
    for (size_t j = 0; j < m.size2(); j++)
      m(i, j) = log(m(i, j));
}

void normalize(matrix_t &m) {
  for (size_t i = 0; i < m.size1(); i++) {
    double z = 0;
    for (size_t j = 0; j < m.size2(); j++)
      z += m(i, j);
    for (size_t j = 0; j < m.size2(); j++)
      m(i, j) /= z;
  }
}

void normalize(vector_t &v) {
  double z = 0;
  for (size_t i = 0; i < v.size(); i++)
    z += v(i);
  for (size_t i = 0; i < v.size(); i++)
    v(i) /= z;
}

matrix_t read_emission(const std::string &path) {
  std::vector<std::vector<double>> d;
  std::string line;
  std::ifstream ifs(path.c_str());
  safeGetline(ifs, line);
  while (not ifs.eof()) {
    size_t idx;
    ifs >> idx;
    double z = 0;
    std::vector<double> v;
    for (size_t i = 0; i < 4; i++) {
      double x;
      ifs >> x;
      v.push_back(x);
      z += x;
    }
    if (fabs(1.0 - z) < 1e-5)
      d.push_back(v);
  }
  matrix_t m(d.size(), d[0].size());
  for (size_t i = 0; i < d.size(); i++)
    for (size_t j = 0; j < d[0].size(); j++)
      m(i, j) = d[i][j];
  return m;
}

confusion_matrix operator+(const confusion_matrix &m, double x) {
  confusion_matrix p(m);
  p.true_positives += x;
  p.false_positives += x;
  p.false_negatives += x;
  p.true_negatives += x;
  return p;
}

string vec2string(const count_vector_t &v) {
  string s = "[" + to_string(v.size()) + "](";
  bool first = true;
  for (auto &x : v) {
    if (first)
      first = false;
    else
      s += ",";
    s += to_string(x);
  }
  s += ")";
  return s;
}
