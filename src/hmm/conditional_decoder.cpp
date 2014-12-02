#include "conditional_decoder.hpp"

using namespace std;

ConditionalDecoder::ConditionalDecoder(const HMM &hmm_) : hmm(hmm_) {
  for(auto &group: hmm.groups) {
    list<matrix_t> matrix_list;
    emission_matrices.push_back(make_pair(group.name, matrix_list));
  }
};

void ConditionalDecoder::decode(std::ostream &os, const Data::Seq &seq) const {
}
