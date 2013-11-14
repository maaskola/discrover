
#include "subhmm.hpp"

using namespace std;

SubHMM::SubHMM(const HMM &hmm, const Training::Range &states) : HMM(hmm, false), n_original_states(n_states), reduce(n_states,-1), lift(states)
{
  sort(lift.begin(), lift.end());
  size_t idx = 0;
  for(auto &state: lift)
    reduce[state] = idx++;

  for(int i = last_state; i >= 0; i--)
    if(find(lift.begin(), lift.end(), i) == lift.end())
      del_column(i);
  initialize_ranges();
  initialize_pred_succ();
  if(verbosity >= Verbosity::debug)
    std::cout << "Constructed SubHMM" << std::endl;
}

Training::Range SubHMM::map_down(const Training::Range &range) const
{
  Training::Range dest;
  for(auto &x: range)
    if(reduce[x] != -1)
      dest.push_back(reduce[x]);
  return(dest);
}

Training::Targets SubHMM::map_down(const Training::Targets &targets) const
{
  Training::Targets dest = {map_down(targets.transition), map_down(targets.emission)};
  return(dest);
}

matrix_t SubHMM::lift_emission(const matrix_t &m) const
{
  matrix_t n = zero_matrix(n_original_states, n_emissions);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      n(lift[i],j) = m(i,j);
  return(n);
}

matrix_t SubHMM::lift_transition(const matrix_t &m) const
{
  matrix_t n = zero_matrix(n_original_states, n_original_states);
  for(size_t i = 0; i < m.size1(); i++)
    for(size_t j = 0; j < m.size2(); j++)
      n(lift[i],lift[j]) = m(i,j);
  return(n);
}

