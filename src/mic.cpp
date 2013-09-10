

#include <set>
#include "mic.hpp"
#include "mic_impl.hpp"

using namespace std;

vector<size_t> optimize_mic(const vector<double> &data, const set<size_t> &fixed_bounds, size_t mic) {
  // std::vector<bool> fixed(data.size()+1, false);
  return(MIC::sample(data, fixed_bounds, mic, 1000000, Verbosity::info));
//  for(auto &i: fixed_bounds)
//    fixed[i] = true;
//  return(vector<size_t>());
}

Data::Collection optimize_mic(const Data::Collection &data, const HMM &hmm, size_t mic) {
  size_t group_idx = 0;
  for(size_t i = 0; i < hmm.get_ngroups(); i++)
    if(hmm.is_motif_group(i)) {
      group_idx = i;
      break;
    }
  set<size_t> fixed_bounds;
  vector<double> occur;
  for(auto &series: data)
    for(auto &set: series) {
      for(auto &seq: set)
        occur.push_back(hmm.posterior_atleast_one(seq, group_idx).posterior);
      fixed_bounds.insert(occur.size());
    }
  // for(size_t i = 0; i < occur.size(); i++)
  //   cout << i << "\t" << occur[i] << endl;
  // for(auto &x: fixed_bounds)
  //   cout << x << endl;
  vector<size_t> bounds = optimize_mic(occur, fixed_bounds, mic);
  auto bounds_iter = bounds.begin();
  size_t cumul = 0;
  Data::Collection mic_data;
  for(auto &series: data) {
    Data::Series ser;
    for(auto &set: series) {
      Data::Set s;
      s.path = set.path;
      s.sha1 = set.sha1;
      s.motifs = set.motifs;
      s.series = set.series;
      for(auto &seq: set) {
        s.sequences.push_back(seq);
        s.set_size++;
        s.seq_size += seq.sequence.size();
        cumul++;
        if(cumul == *bounds_iter) {
          ser.sets.push_back(s);
          ser.set_size += s.set_size;
          ser.seq_size += s.seq_size;
          s = Data::Set();
          s.path = set.path;
          s.sha1 = set.sha1;
          s.motifs = set.motifs;
          s.series = set.series;
          bounds_iter++;
        }
      }
    }
    mic_data.series.push_back(ser);
    mic_data.set_size += ser.set_size;
    mic_data.seq_size += ser.seq_size;
  }

  return(mic_data);
}


