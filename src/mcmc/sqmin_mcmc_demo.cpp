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
 *       Filename:  sqmin_mcmc_demo.hpp
 *
 *    Description:  Demonstration how to use the MCMC code
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#include <ctime>
#include <iostream>
#include "montecarlo.hpp"

namespace MCMC {
  template <>
    class Evaluator<double> {
      public:
        double evaluate(double i) const { return -i*i; };
    };
  template <>
    class Generator<double> {
      public:
        double generate(double i) const { return i+2*(rand() * 1.0 / RAND_MAX -0.5); };
    };
}

using namespace std;

int main(int argc, char** argv) {
  srand(time(0));

  size_t n_iter = 10;
  const size_t n_chains = 10;
  double temperature = 1e-3;

  if(argc > 1)
    n_iter = atoi(argv[1]);

  cout << "Doing " << n_iter << " iterations." << endl;

  MCMC::Evaluator<double> eval;
  MCMC::Generator<double> gen;
  MCMC::MonteCarlo<double> mcmc(gen, eval, Verbosity::info);
  vector<double> temperatures;
  vector<double> init;
  for(size_t i = 0; i < n_chains; i++) {
    init.push_back(100 + i); // + 1.0 * rand() / RAND_MAX);
    temperatures.push_back(temperature);
    temperature /= 2;
  }
  cout << "Initial values:";
  for(auto &x: init)
    cout << "\t" << x;
  cout << endl;
  auto results = mcmc.parallel_tempering(temperatures, init, n_iter);
  cout << "Final values:";
  for(auto &x: results)
    cout << "\t" << x.back().first;
  cout << endl;
  return(EXIT_SUCCESS);
}

