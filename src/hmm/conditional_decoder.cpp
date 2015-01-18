#include "conditional_decoder.hpp"

using namespace std;

ConditionalDecoder::ConditionalDecoder(const HMM &hmm_) : hmm(hmm_) {
  for (auto &group : hmm.groups)
    if (group.kind == Group::Kind::Motif) {
      vector<matrix_t> matrices;
      // TODO: handle indels in the motifs
      size_t n = group.states.size();
      matrix_t matrix(n, HMM::n_emissions);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < HMM::n_emissions; ++j)
          matrix(i, j) = hmm.emission(group.states[0] + i, j);
      matrices.push_back(matrix);

      emission_matrices.push_back(make_pair(group.name, matrices));
    }

  if (false) {
    cerr << "Conditional Decoder" << endl;
    for (auto &x : emission_matrices)
      for (auto &y : x.second)
        cerr << x.first << endl << y << endl;
  }
};

void ConditionalDecoder::decode(std::ostream &os, const Data::Seq &seq) const {
  const size_t n = seq.sequence.size();
  for (auto &group : emission_matrices)
    for (auto &matrix : group.second) {
      // TODO: handle indels in the motifs
      // -> number the emission matrix variants
      os << "Conditional (" << group.first << ")";
      const size_t w = matrix.size1();
      for (size_t i = 0; i < n - w + 1; ++i) {
        double p = 1;
        for (size_t j = 0; j < w; ++j)
          p *= matrix(j, seq.isequence(i + j));
        os << " " << p;
      }
      for (size_t i = 1; i < w; ++i)
        os << " " << 0;
      os << endl;
    }
}
