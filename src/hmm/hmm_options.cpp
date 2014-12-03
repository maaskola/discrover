#include "hmm_options.hpp"
#include "trainingmode.hpp"

using namespace std;

namespace Options {

  string compression2ending(Compression compression)
  {
    switch(compression) {
      case Compression::none:
        return("");
      case Compression::gzip:
        return(".gz");
      case Compression::bzip2:
        return(".bz2");
    }
    return("");
  }

  string compression2string(Compression compression)
  {
    switch(compression) {
      case Compression::none:
        return("none");
      case Compression::gzip:
        return("gzip");
      case Compression::bzip2:
        return("bzip2");
    }
    return("");
  }

  istream& operator>>(istream &in, Compression &compression)
  {
    string token;
    in >> token;
    if(token == "none")
      compression = Compression::none;
    else if(token == "gzip" or token == "gz")
      compression = Compression::gzip;
    else if(token == "bzip2" or token == "bz2")
      compression = Compression::bzip2;
    else {
      cout << "Error: can not parse compression type '" << token << "'." << endl;
      exit(-1);
    }
    return in;
  }

  istream &operator>>(istream &is, MultiMotif::Relearning &relearning)
  {
    string word;
    is >> word;
    if(word == "none")
      relearning = MultiMotif::Relearning::None;
    else if(word == "reest")
      relearning = MultiMotif::Relearning::Reestimation;
    else if(word == "full")
      relearning = MultiMotif::Relearning::Full;
    else {
      cout << "Error: can not parse relearning mode '" << word << "'." << endl;
      exit(-1);
    }
    return(is);
  }

  ostream& operator<<(ostream &out, const Compression &compression)
  {
    out << compression2string(compression);
    return out;
  }

  ostream &operator<<(ostream &os, const MultiMotif::Relearning &relearning)
  {
    switch(relearning) {
      case MultiMotif::Relearning::None:
        os << "none";
        break;
      case MultiMotif::Relearning::Reestimation:
        os << "reest";
        break;
      case MultiMotif::Relearning::Full:
        os << "full";
        break;
    }
    return(os);
  }

  ostream &operator<<(ostream &os, const Verbosity &verbosity)
  {
    switch(verbosity) {
      case Verbosity::nothing:
        os << "nothing" << endl;
        break;
      case Verbosity::error:
        os << "error" << endl;
        break;
      case Verbosity::info:
        os << "info" << endl;
        break;
      case Verbosity::verbose:
        os << "verbose" << endl;
        break;
      case Verbosity::debug:
        os << "debug" << endl;
        break;
      case Verbosity::everything:
        os << "everything" << endl;
        break;
    }
    return(os);
  }

  template <typename X> ostream& operator<<(ostream &os, const vector<X> &v)
  {
    bool first = true;
    for(auto &x: v) {
      if(first)
        first = false;
      else
        os << " ";
      os << x;
    }
    return(os);
  }

  ostream &operator<<(ostream &os, const LineSearch &options)
  {
    os << "Line search options:" << endl
      << "mu = " << options.mu << endl
      << "eta = " << options.eta << endl
      << "delta = " << options.delta << endl
      << "max_steps = " << options.max_steps << endl;
    return(os);
  }

  ostream &operator<<(ostream &os, const Termination &options)
  {
    os << "Termination options:" << endl
      << "max_iter = " << options.max_iter << endl
      << "past = " << options.past << endl
      << "gamma_tolerance = " << options.gamma_tolerance << endl
      << "delta_tolerance = " << options.delta_tolerance << endl
      << "epsilon_tolerance = " << options.epsilon_tolerance << endl
      << "absolute_improvement = " << options.absolute_improvement << endl;
    return(os);
  }

  ostream &operator<<(ostream &os, const Sampling &options)
  {
    os << "Sampling options:" << endl
      << "do_sampling = " << options.do_sampling << endl
      << "min_size = " << options.min_size << endl
      << "max_size = " << options.max_size << endl
      << "temperature = " << options.temperature << endl
      << "anneal_factor = " << options.anneal_factor << endl
      << "n_indels = " << options.n_indels << endl
      << "n_shift = " << options.n_shift << endl
      << "n_parallel = " << options.n_parallel << endl;
    return(os);
  }

  ostream &operator<<(ostream &os, const ExecutionInformation &exec_info)
  {
    os << "ExecutionInformation:" << endl << "STUB!" << endl;
    return(os);
  }

  ostream &operator<<(ostream &os, const Evaluation &eval_info)
  {
    os << "Evaluate, perform RIC analysis = " << eval_info.perform_ric << endl;
    os << "Evaluate, skip Viterbi path = " << eval_info.skip_viterbi_path << endl;
    os << "Evaluate, skip occurrence table = " << eval_info.skip_occurrence_table << endl;
    os << "Evaluate, skip summary = " << eval_info.skip_summary << endl;
    os << "Evaluate, print_posterior = " << eval_info.print_posterior << endl;
    return(os);
  }

  ostream &operator<<(ostream &os, const HMM &options)
  {
    os << "HMM options:" << endl
      << "paths = " << options.paths << endl
      << "motif_specifications = " << options.motif_specifications << endl
      << "load_paths = " << options.load_paths << endl
      << "label = " << options.label << endl
      << "dont_save_shuffle_sequences = " << options.dont_save_shuffle_sequences << endl
      << "seeds = " << options.seeds << endl
      // << "plasma_options = " << options.plasma_options << endl // TODO: implement
      << "evaluation_options = " << options.evaluate << endl
      << "n_threads = " << options.n_threads << endl
      << "n_seq = " << options.n_seq << endl
      << "alpha = " << options.alpha << endl
      << "contingency_pseudo_count = " << options.contingency_pseudo_count << endl
      << "emission_pseudo_count = " << options.emission_pseudo_count << endl
      << "transition_pseudo_count = " << options.transition_pseudo_count << endl
      << "n_simulations = " << options.n_simulations << endl
      << "long_names = " << options.long_names << endl
      << "class_model = " << options.class_model << endl
      << "revcomp = " << options.revcomp << endl
      << "weighting = " << options.weighting << endl
      << "accept_multiple = " << options.multi_motif.accept_multiple << endl
      << "relearning = " << options.multi_motif.relearning << endl
      << "residual_ratio = " << options.multi_motif.residual_ratio << endl
      << "output_compression = " << options.output_compression << endl
      << "extend= " << options.extend << endl
      << "left_padding = " << options.left_padding << endl
      << "right_padding = " << options.right_padding << endl
      << "timing_information = " << options.timing_information << endl
      << "cross_validation_iterations = " << options.cross_validation_iterations << endl
      << "cross_validation_freq = " << options.cross_validation_freq << endl
      << "store_intermediate = " << options.store_intermediate << endl
      << "wiggle = " << options.wiggle << endl
      << "line_search = " << options.line_search << endl
      << "random_salt = " << options.random_salt << endl
      << "dont_learn_class_prior = " << options.dont_learn_class_prior << endl
      << "dont_learn_conditional_motif_prior = " << options.dont_learn_conditional_motif_prior << endl
      << "class_prior = " << options.class_prior << endl
      << "conditional_motif_prior1 = " << options.conditional_motif_prior1 << endl
      << "conditional_motif_prior2 = " << options.conditional_motif_prior2 << endl
      << "bg_learning = " << options.bg_learning << endl
      // << "objectives = " << options.objectives << endl // TODO: implement
      << "termination = " << options.termination << endl
      << "limit_logp = " << options.limit_logp << endl
      << "miseeding = " << options.use_mi_to_seed<< endl
      << "sampling = " << options.sampling << endl
      << "lambda = " << options.lambda << endl
      << "emission_matrix_paths = " << options.emission_matrix_paths << endl
      << "verbosity = " << options.verbosity << endl
      << "exec_info = " << options.exec_info << endl;
    return(os);
  }

};
