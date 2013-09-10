#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>

double rgamma( double mean, double variance, boost::mt19937& rng ) {
  const double shape = ( mean*mean )/variance;
  double scale = variance/mean;

  boost::gamma_distribution<> gd( shape );
  boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );

  return scale*var_gamma();
}

std::vector<double> rdirichlet(std::vector<double> alpha, boost::mt19937& rng)
{
  size_t n = alpha.size();
  std::vector<double> v(n);
  double z = 0;
  for(size_t i = 0; i < n; i++)
    z += v[i] = rgamma(alpha[i], 1.0, rng);
  for(size_t i = 0; i < n; i++)
    v[i] /= z;
  return(v);
}

std::vector<double> runiform(size_t n, boost::mt19937& rng)
{
  std::vector<double> v(n);
  double z = 0;
  for(size_t i = 0; i < n; i++)
    z += v[i] = rgamma(1.0, 1.0, rng);
  for(size_t i = 0; i < n; i++)
    v[i] /= z;
  return(v);
}

/*
template<class C> std::ostream &operator<<(std::ostream &os, const std::vector<C> &v)
{
  for(size_t i = 0; i < v.size(); i++)
    os << (i!=0 ? " " : "") << v[i];
  return(os);
}
*/

/*
//#include <iostream>
//#include <string>

using namespace std;

int main(int argc, char **argv)
{
  size_t n = atoi(argv[1]);
  size_t k = atoi(argv[2]);
  boost::mt19937 rng;
  for(size_t i = 0; i < n; i++)
    cout << runiform(k, rng) << endl;
  return(EXIT_SUCCESS);
}
*/

