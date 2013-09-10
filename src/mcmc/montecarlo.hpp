
#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP


#include <cstdlib>
#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include "../verbosity.hpp"

namespace MCMC {

  double boltzdist(double dG, double T){
    return exp(- dG / T);
  };

  template <class T>
    class Evaluator {
      public:
        double evaluate(T& i) const;
    };

  template <>
    class Evaluator<double> {
      public:
        double evaluate(double i) const { return i*i; };
    };

  template <class T>
    class Generator {
      public:
        T generate(const T& i) const;
    };

  template <>
    class Generator<double> {
      public:
        double generate(double i) const { return i+2*(rand() * 1.0 / RAND_MAX -0.5); };
    };


  template <class T>
    class MonteCarlo {
      typedef std::pair<T,double> E;
      Verbosity verbosity;
      public:
        MonteCarlo(Verbosity ver) :
          verbosity(ver),
          generator(Generator<T>()),
          evaluator(Evaluator<T>()) { };
        MonteCarlo(const Generator<T> &gen, const Evaluator<T> &eval, Verbosity ver):
          verbosity(ver),
          generator(gen),
          evaluator(eval) { };
        ~MonteCarlo(){};

        Generator<T> generator;
        Evaluator<T> evaluator;

      private:
        bool GibbsStep(double temp, T &state, double &G) const {
          T nextstate = generator.generate(state);
          double nextG = evaluator.evaluate(nextstate);
          double dG = nextG - G;
          double r = 1.0 * rand() / RAND_MAX;
          double p = std::min<double>(1.0, boltzdist(-dG,temp));
          if(verbosity >= Verbosity::verbose)
            std::cerr << "T = " << temp << " next state = " << nextstate << std::endl
              << "nextG = " << nextG << " G = " << G << " dG = " << dG << std::endl
              << " p = " << p << " r = " << r << std::endl;
          if(isnan(nextG) == 0 and (dG > 0 or r <= p)) {
            if(verbosity >= Verbosity::verbose)
              std::cerr << "Accepted!" << std::endl;
            state = nextstate;
            G = nextG;
            return(true);
          } else {
            if(verbosity >= Verbosity::verbose)
              std::cerr << "Rejected!" << std::endl;
            return(false);
          }
       }
 

        bool swap(double temp1, double temp2, T &state1, T &state2, double &G1, double &G2) const {
          double r = 1.0 * rand() / RAND_MAX;
          double p = std::min<double>(1.0,exp(-(G1 / temp1 + G2 / temp2 - G1 / temp2 - G2 / temp1)));
          if(verbosity >= Verbosity::verbose)
            std::cerr << "T1 = " << temp1 << " T2 " << temp2 << " G1 = " << G1 << " G2 = " << G2 << std::endl
              << "r = " << r << " p = " << p << std::endl;
          if(r <= p) {
            if(verbosity >= Verbosity::verbose)
              std::cerr << "Swap!" << std::endl;
            std::swap<T>(state1, state2);
            std::swap<double>(G1, G2);
            return true;
          } else {
            if(verbosity >= Verbosity::verbose)
              std::cerr << "Swap rejected!" << std::endl;
            return false;
          }
        }


      public:
        std::list<E> run(double temp, double anneal, const T &init, size_t steps) const {
          T state = T(init);
          double G = evaluator.evaluate(state);
          std::list<E> trajectory;
          trajectory.push_back(E(state,G));
          for(size_t i = 0; i < steps; i++) {
            if(GibbsStep(temp, state, G))
              trajectory.push_back(E(state,G));
            temp *= anneal;
          }
          return trajectory;
        };


        std::vector<std::list<E>> parallel_tempering(const std::vector<double> &temp, const std::vector<T> &init, size_t steps) const {
          size_t n = temp.size();
          std::vector<T> state = init;

          std::vector<double> G;
          for(auto s : state)
            G.push_back(evaluator.evaluate(s));

          std::vector<std::list<E>> trajectory(temp.size());
          for(size_t t = 0; t < temp.size(); t++)
            trajectory[t].push_back(E(state[t],G[t]));

          for(size_t i = 0; i < steps; i++) {
            if(verbosity >= Verbosity::info)
              std::cerr << "Iteration " << i << " of " << steps << std::endl;
            for(size_t t = 0; t < n; t++)
              if(GibbsStep(temp[t], state[t], G[t])) // TODO: if one wants to determine means one should respect the failed changes, and input once more the original state to the trajectory.
                trajectory[t].push_back(E(state[t],G[t]));

            if(verbosity >= Verbosity::info) {
              std::cout << "Scores =";
              for(size_t t = 0; t < n; t++)
                std::cout << " " << G[t];
              std::cout << std::endl;
            }

            if(temp.size() > 1) {
              size_t r = rand() % (n - 1);
              if(verbosity >= Verbosity::verbose)
                std::cerr << "Testing swap of " << r << " and " << r+1 << std::endl;
              if(swap(temp[r], temp[r+1], state[r], state[r+1], G[r], G[r+1])) {
                trajectory[r].push_back(E(state[r],G[r]));
                trajectory[r+1].push_back(E(state[r+1],G[r+1]));
                if(verbosity >= Verbosity::info)
                  std::cerr << "Swapping chains " << r << " and " << r+1 << "." << std::endl;
              }
            }
          }
          return trajectory;
        };
    };
}

#endif

