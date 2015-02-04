
#include <iomanip>
#include <fstream>
#include <vector>
#include <numeric>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include "../aux.hpp"
#include "report.hpp"
#include "conditional_decoder.hpp"
#include "../timer.hpp"
#include "../format_constants.hpp"
#include "../plasma/plasma.hpp"

#if LIBR_FOUND
  #define MATHLIB_STANDALONE
  #include <Rmath.h>
#else
  #include "../stats/pgamma.hpp"
#endif

using namespace std;

#if CAIRO_FOUND
#include "../logo/logo.hpp"

vector<string> Evaluator::generate_logos(const string &path_stem,
                                         const Options::HMM &options) const {
  vector<string> paths;
  size_t motif_idx = 0;
  for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
    if (hmm.is_motif_group(group_idx)) {
      const string nucls = "acgt";
      Logo::matrix_t matrix;
      for (auto state : hmm.groups[group_idx].states) {
        Logo::column_t col(4, 0);
        for (size_t i = 0; i < nucls.size(); i++)
          col[i] = hmm.emission(state, i);
        matrix.push_back(col);
      }
      for (auto path : Logo::draw_logo(matrix, path_stem, options.logo))
        paths.push_back(path);
      motif_idx++;
    }
  return paths;
}
#endif

string contrast_name_tag(const Data::Contrast &contrast) {
  string tag = contrast.name + ":";
  bool first = true;
  for (auto &dataset : contrast) {
    tag += (first ? "" : ",") + dataset.name();
    first = false;
  }
  return tag;
}

void print_table(ostream &ofs, const matrix_t m, const Data::Contrast &contrast,
                 size_t width, size_t prec) {
  const ios::fmtflags flags(ofs.flags());
  size_t w = 0;
  for (size_t i = 0; i < m.size1(); i++)
    w = max(w, contrast.sets[i].name().size());

  ofs << left << setw(w) << "" << right << setw(width) << "Present" << right
      << setw(width) << "Absent";
  if (m.size2() == 2)
    ofs << right << setw(width) << "Percent";
  ofs << endl;

  for (size_t i = 0; i < m.size1(); i++) {
    ofs << left << setw(w) << contrast.sets[i].name();
    for (size_t j = 0; j < m.size2(); j++)
      ofs << right << fixed << setprecision(prec) << setw(width) << m(i, j);
    if (m.size2() == 2)
      ofs << right << fixed << setprecision(prec) << setw(width)
          << (m(i, 0) / (m(i, 0) + m(i, 1)) * 100);
    ofs << left << endl;
  }
  ofs.flags(flags);
}

template <typename X>
void print(ostream &ofs, const string &tag, const string &what, X x,
           const string &unit = "") {
  ofs << tag << what << " = " << x;
  if (unit != "")
    ofs << " " << unit;
  ofs << endl;
}

void count_report(ostream &ofs, const matrix_t counts, size_t motif_len,
                  const Data::Contrast &contrast, double pseudo_count,
                  bool limit_logp, const string &name, const string &tag) {
  size_t degrees_freedom = (counts.size1() - 1) * (counts.size2() - 1);
  print_table(ofs, counts, contrast, 10, 2);

  double mutualinf
      = calc_mutual_information(counts, pseudo_count, true, false, false);
  double exp_mutualinf
      = calc_mutual_information(counts, pseudo_count, true, true, false);
  double var_mutualinf
      = calc_mutual_information(counts, pseudo_count, true, true, true);
  double sd_mutualinf = sqrt(var_mutualinf);
  double zscore_mutualinf = exp_mutualinf / sd_mutualinf;
  double g = calc_g_test(counts, pseudo_count);
  // TODO: reactivate
  // double correct_class = hmm.correct_classification(contrast);
  double log_p_g = pchisq(g, degrees_freedom, false, true);
  double cor_log_p_g_stringent = log(149) * motif_len + log_p_g;
  if (limit_logp)
    cor_log_p_g_stringent = min<double>(0, cor_log_p_g_stringent);

  print(ofs, tag, "Discriminatory mutual information", mutualinf,
        "bit per sequence");
  print(ofs, tag, "Expected discriminatory mutual information", exp_mutualinf,
        "bit per sequence");
  print(ofs, tag, "Variance of discriminatory mutual information",
        var_mutualinf, "bit per sequence");
  print(ofs, tag, "Std. dev. of discriminatory mutual information",
        sd_mutualinf, "bit per sequence");
  print(ofs, tag, "Z-score of discriminatory mutual information",
        zscore_mutualinf, "bit per sequence");
  print(ofs, tag, "G-test", g);
  print(ofs, tag, "Log-P(Chi-Square(G-Test))", log_p_g);
  print(ofs, tag, "Bonferroni corrected log-P(Chi-Square(G-Test))",
        cor_log_p_g_stringent);
}

Evaluator::Evaluator(const HMM &hmm_) : hmm(hmm_){};

void Evaluator::eval_contrast(ostream &ofs, const Data::Contrast &contrast,
                              bool limit_logp, const string &tag) const {
  // double mi = hmm.mutual_information(contrast, contrast);
  // ofs << "Summed discriminatory mutual information = " << mi << " bit per
  // sequence" << endl;

  for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
    if (hmm.is_motif_group(group_idx)) {
      size_t motif_len = hmm.get_motif_len(group_idx);
      bitmask_t present_mask = 1 << group_idx;
      // TODO EVALUATION
      vector_t v = hmm.posterior_atleast_one(contrast, present_mask);
      matrix_t counts(v.size(), 2);
      for (size_t i = 0; i < v.size(); i++) {
        counts(i, 0) = v(i);
        counts(i, 1) = contrast.sets[i].set_size - v(i);
      }
      // TODO EVALUATION
      // TODO compute from count table
      double mcc = hmm.matthews_correlation_coefficient(contrast, present_mask);
      // double llr = calc_log_likelihood_ratio(counts, hmm.pseudo_count);
      // double dips_tscore = hmm.dips_tscore(contrast, feature);
      // TODO EVALUATION
      double dips_sitescore = hmm.dips_sitescore(contrast, present_mask);
      // TODO: reactivate
      // double correct_class = hmm.correct_classification(contrast);

      string name = hmm.get_group_name(group_idx);
      string consensus = hmm.get_group_consensus(group_idx);

      ofs << tag << "Expected count of sequences with at least one occurrence "
                    "of motif '" + name + ":" + consensus << "'" << endl;
      count_report(ofs, counts, motif_len, contrast, hmm.get_pseudo_count(),
                   limit_logp, name, tag);
      print(ofs, tag, "Matthews correlation coefficient", mcc);
      // print(ofs, tag, "DIPS t-score", dips_tscore * 100, "%");
      print(ofs, tag, "DIPS site-score", dips_sitescore * 100, "%");
      // TODO: reactivate
      // print(ofs, tag, "log P correct classification", correct_class);
      // TODO print a table of the likelihoods

      // TODO EVALUATION
      // TODO compute from count table
      print(ofs, tag, "Log likelihood difference",
            hmm.log_likelihood_difference(contrast, present_mask));
    }
}

template <class X, class Y>
double cor_pearson(const vector<X> &x, const vector<Y> &y) {
  size_t n = x.size();
  double zx = accumulate(x.begin(), x.end(), 0);
  double zy = accumulate(y.begin(), y.end(), 0);
  zx /= n;
  zy /= n;

  double num = 0;
  double denom_x = 0, denom_y = 0;
  for (size_t i = 0; i < n; i++) {
    double a = x[i] - zx;
    double b = y[i] - zy;
    num += a * b;
    denom_x += a * a;
    denom_y += b * b;
  }
  if (num == 0)
    return 0;
  double r = num / sqrt(denom_x) / sqrt(denom_y);
  return r;
}

template <class X>
vector<double> tied_ranking(const vector<X> &x) {
  const size_t n = x.size();
  vector<double> r(n);
  map<X, size_t> counts;
  for (auto &y : x)
    counts[y]++;
  size_t cumul = 0;
  for (auto &y : counts) {
    double first = cumul;
    cumul += y.second;
    double ave = (first + cumul - 1) / 2;
    y.second = ave;
  }
  for (size_t i = 0; i < n; i++)
    r[i] = counts[x[i]];
  return r;
}

template <class X, class Y>
double cor_spearman(const vector<X> &x, const vector<Y> &y) {
  return cor_pearson(tied_ranking(x), tied_ranking(y));
}

template <class X>
double cor_rank(const vector<X> &x) {
  vector<X> y(x.size());
  iota(y.begin(), y.end(), 0);
  return cor_pearson(tied_ranking(x), y);
}

double cor_fisher_z(double x, size_t n) {
  double z = sqrt((n - 3) / 1.06) * atanh(x);
  return z;
}

double cor_student_t(double x, size_t n) {
  double t = x * sqrt((n - 2) / (1 - x * x));
  return t;
}

string stars(double x) {
  if (x < 0.001) return "***";
  if (x < 0.01)  return "**";
  if (x < 0.05)  return "*";
  return "";
}

const boost::math::normal_distribution<double> standard_normal_distribution;

template <class T>
void correlation_report(const vector<T> &x, ostream &out, size_t width = 12,
                        size_t prec = 5) {
  const ios::fmtflags flags(out.flags());
  double rho = cor_rank(x);
  size_t n = x.size();
  double z = cor_fisher_z(rho, n);
  double t = cor_student_t(rho, n);
  double p_norm = 0, p_t = 0;
  {
    using namespace boost::math;
    const students_t_distribution<> students_t_dist(n - 2);
    if (rho < 0) {
      p_norm = cdf(standard_normal_distribution, z);
      p_t = cdf(students_t_dist, t);
    } else {
      p_norm = cdf(complement(standard_normal_distribution, z));
      p_t = cdf(complement(students_t_dist, t));
    }
  }
  out << setw(width) << fixed << right << setprecision(prec) << rho
      << setw(width) << fixed << right << setprecision(2) << z << setw(width)
      << fixed << right << setprecision(2) << log(p_norm) << setw(4) << right
      << stars(p_norm) << setw(width) << fixed << right << setprecision(2) << t
      << setw(width) << fixed << right << setprecision(2) << log(p_t) << setw(4)
      << right << stars(p_t) << endl;
  out.flags(flags);
}

void Evaluator::print_posterior(ostream &os, const vector_t &scale,
                                const matrix_t &f, const matrix_t &b) const {
  const size_t n = scale.size() - 2;
  vector_t posterior(n);
  for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
    if (hmm.is_motif_group(group_idx)) {
      size_t k = *begin(hmm.groups[group_idx].states);
      os << "Posterior (" << hmm.get_group_name(group_idx) << ")";
      for (size_t i = 1; i <= n; ++i)
        os << " " << f(i, k) * b(i, k) * scale(i);
      os << endl;
    }
}

Evaluator::ResultsCounts Evaluator::evaluate_dataset(
    const Data::Set &dataset, ostream &out, ostream &v_out, ostream &occ_out,
    ostream &bed_out, const Options::HMM &options) const {
  const size_t width = 12;
  const size_t prec = 5;
  Timer timer;

  ConditionalDecoder conditional_decoder(hmm);

  map<size_t, double> n_sites;
  map<size_t, double> n_motifs;
  map<size_t, size_t> n_viterbi_sites;
  map<size_t, size_t> n_viterbi_motifs;

  if (not options.evaluate.skip_summary) {
    // TODO EVALUATION
    double log_likelihood = hmm.log_likelihood(dataset);
    // out << endl << "Summary of " << dataset.name() << endl;
    // out << "Total sequences = " << dataset.sequences.size() << endl;
    out << "Log-likelihood of " << dataset.name() << " = " << log_likelihood
        << endl;
    // out << "Akaike information criterion AIC = " << 2 * hmm.n_parameters() -
    // 2 * log_likelihood << endl;
    //      out << "Akaike information criterion AIC = " << 2 *
    //      hmm.non_zero_parameters(hmm.gen_training_targets(options)) - 2 *
    //      log_likelihood << endl;
  }

  if (not options.evaluate.skip_viterbi_path)
    v_out << "# " << dataset.name() << " details following" << endl;

  const size_t n_groups = hmm.get_ngroups();
  const size_t n = dataset.sequences.size();
  size_t number_motifs = 0;
  for (size_t group_idx = 0; group_idx < n_groups; group_idx++)
    if (hmm.is_motif_group(group_idx))
      number_motifs++;

  vector<vector<double>> atl_counts(number_motifs), exp_counts(number_motifs),
      vit_counts(number_motifs);
  for (size_t group_idx = 0; group_idx < number_motifs; group_idx++) {
    atl_counts[group_idx] = vector<double>(n);
    exp_counts[group_idx] = vector<double>(n);
    vit_counts[group_idx] = vector<double>(n);
  }

  if (options.evaluate.perform_ric)
    for (size_t group_idx = 0; group_idx < n_groups; group_idx++)
      if (hmm.is_motif_group(group_idx)) {
        bitmask_t present_mask = 1 << group_idx;
        // TODO EVALUATION
        double ric = hmm.rank_information(dataset, present_mask);
        out << "RIC = " << ric << endl;
      }

  for (size_t i = 0; i < dataset.sequences.size(); i++) {
    HMM::StatePath viterbi_path;
    // TODO EVALUATION
    double lp = hmm.viterbi(dataset.sequences[i], viterbi_path);

    if (not(options.evaluate.skip_viterbi_path
            and options.evaluate.skip_summary
            and options.evaluate.skip_bed)) {
      stringstream viterbi_str, exp_str, atl_str;

      bool first = true;
      size_t motif_idx = 0;
      for (size_t group_idx = 0; group_idx < n_groups; group_idx++)
        if (hmm.is_motif_group(group_idx)) {
          if (first)
            first = false;
          else {
            viterbi_str << "/";
            exp_str << "/";
            atl_str << "/";
          }

          bitmask_t present_mask = 1 << group_idx;
          // TODO EVALUATION
          double atl = hmm.posterior_atleast_one(dataset.sequences[i],
                                                 present_mask).posterior;
          // TODO EVALUATION
          double expected
              = hmm.expected_posterior(dataset.sequences[i], present_mask);
          size_t n_viterbi = hmm.count_motif(viterbi_path, group_idx);
          atl_counts[motif_idx][i] = atl;
          exp_counts[motif_idx][i] = expected;
          vit_counts[motif_idx][i] = n_viterbi;

          n_sites[group_idx] += atl;
          n_motifs[group_idx] += expected;
          n_viterbi_sites[group_idx] += (n_viterbi > 0 ? 1 : 0);
          n_viterbi_motifs[group_idx] += n_viterbi;

          if (not options.evaluate.skip_viterbi_path) {
            viterbi_str << n_viterbi;
            exp_str << expected;
            atl_str << atl;
          }

          motif_idx++;
        }
      if (not options.evaluate.skip_viterbi_path) {
        v_out << ">" << dataset.sequences[i].definition << endl;
        v_out << "V-sites = " << viterbi_str.str()
              << " E-sites = " << exp_str.str()
              << " P(#sites>=1) = " << atl_str.str()
              << " Viterbi log-p = " << lp << endl;
        v_out << dataset.sequences[i].sequence << endl;
        v_out << hmm.path2string_group(viterbi_path) << endl;
      }

      if (options.evaluate.print_posterior) {
        vector_t scale;
        auto f = hmm.compute_forward_scaled(dataset.sequences[i], scale);
        auto b = hmm.compute_backward_prescaled(dataset.sequences[i], scale);
        print_posterior(v_out, scale, f, b);
      }
      if (options.evaluate.conditional_motif_probability)
        conditional_decoder.decode(v_out, dataset.sequences[i]);
    }

    if (not options.evaluate.skip_bed)
      hmm.print_occurrence_table(dataset.name(), dataset.sequences[i],
                                 viterbi_path, bed_out, true);
    if (not options.evaluate.skip_occurrence_table)
      hmm.print_occurrence_table(dataset.name(), dataset.sequences[i],
                                 viterbi_path, occ_out, false);
  }

  if (not options.evaluate.skip_summary) {
    const size_t col_width = 17;
    out << setw(col_width) << left << "Rank analysis";
    out << setw(width) << right << "Rho" << setw(width) << right << "Z"
        << setw(width) << right << "log P(Z)" << setw(4) << right << ""
        << setw(width) << right << "t" << setw(width) << right << "log P(t)"
        << setw(4) << right << "" << endl;
    // TODO: write out the motif name
    for (size_t group_idx = 0; group_idx < number_motifs; group_idx++) {
      out << setw(col_width) << left << "Expected sites";
      correlation_report(atl_counts[group_idx], out, width, prec);
      out << setw(col_width) << left << "Expected motifs";
      correlation_report(exp_counts[group_idx], out, width, prec);
      vector<size_t> vit(n);
      for (size_t i = 0; i < n; i++)
        vit[i] = vit_counts[group_idx][i] > 0;
      out << setw(col_width) << left << "Viterbi sites";
      correlation_report(vit, out, width, prec);
      out << setw(col_width) << left << "Viterbi motifs";
      correlation_report(vit_counts[group_idx], out, width, prec);
    }
  }

  double time = timer.tock();
  if (options.timing_information)
    cerr << "Evaluation for " + dataset.name() + ": "
            + time_to_pretty_string(time) << endl;
  ResultsCounts results
      = {n_sites, n_motifs, n_viterbi_sites, n_viterbi_motifs};
  return results;
}

Evaluator::Result Evaluator::report(const Data::Collection &collection,
                                    const string &tag,
                                    const Training::Tasks &tasks,
                                    const Options::HMM &options) const {
  Result result;
  // TODO see that this does not invalidate previously learned parameters for
  // MMIE!
  // for(auto &contrast: collection)
  //   for(auto &dataset: contrast)
  //     register_dataset(dataset,
  //     (1.0*dataset.set_size)/training_data.set_size,
  //     options.conditional_motif_prior1, options.conditional_motif_prior2);
  /*  TODO: re-enable class-based models
      if(options.adapt_classes) {
      if(options.verbosity >= Verbosity::info)
      cout << "Registering collection sets for class based HMMs." << endl;

      for(auto &contrast: collection)
      for(auto &dataset: contrast)
      hmm.register_class_hmm(dataset);
      }
      Training::Targets generative_targets =
     hmm.gen_secondary_training_targets(options.training_type);
      Training::Task generative_task;
      generative_task.method = Training::Method::reestimation;
      generative_task.measure = Measure::likelihood;
      generative_task.targets = generative_targets;
      hmm.reestimation(dataset, generative_task, options);
      */

  ios_base::openmode flags = ios_base::out;
  if (options.output_compression != Options::Compression::none)
    flags |= ios_base::binary;

  string file_tag = "";
  if (tag != "")
    file_tag = tag + ".";
  result.files.summary = options.label + file_tag + ".summary";
  result.files.table = options.label + file_tag + ".table"
                       + compression2ending(options.output_compression);
  result.files.viterbi = options.label + file_tag + ".viterbi"
                         + compression2ending(options.output_compression);
  result.files.bed = options.label + file_tag + ".bed"
                     + compression2ending(options.output_compression);

#if CAIRO_FOUND
  result.files.logos = generate_logos(options.label + file_tag, options);
#endif

  ofstream summary_out, occurrence_file, viterbi_file, bed_file;
  summary_out.open(result.files.summary.c_str());

  if (not options.evaluate.skip_summary) {
    const size_t base_col_width = 10;
    size_t col0_w = base_col_width, col1_w = base_col_width,
           col2_w = base_col_width;
    for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
      // do not output information about non-motif groups
      if (hmm.is_motif_group(group_idx)) {
        string consensus = hmm.get_group_consensus(group_idx);
        string name = hmm.get_group_name(group_idx);
        col0_w = max(col0_w, name.size());
        col1_w = max(col1_w, consensus.size());
      }
    col0_w += 2;
    col1_w += 2;
    summary_out << "Motif summary" << endl;
    // TODO: print out PWM and affinity matrix
    summary_out << left << setw(col0_w) << "Motif name" << left << setw(col1_w)
                << "Consensus" << right << setw(col2_w) << "IC [bit]" << endl;
    for (size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
      // do not output information about non-motif groups
      if (hmm.is_motif_group(group_idx)) {
        string consensus = hmm.get_group_consensus(group_idx);
        string name = hmm.get_group_name(group_idx);
        double ic = hmm.information_content(group_idx);
        summary_out << left << setw(col0_w) << name << left << setw(col1_w)
                    << consensus << right << setw(col2_w) << ic << endl;
      }
    summary_out << endl;

    Timer eval_timer;
    for (auto &contrast : collection) {
      summary_out << endl << "Discriminative statistics for contrast '"
                  << contrast_name_tag(contrast) << "'" << endl;
      eval_contrast(summary_out, contrast, options.limit_logp, tag);
    }
    double time = eval_timer.tock();

    if (options.timing_information)
      cerr << "Evaluation of contrast for " + (tag == "" ? "" : tag + " ")
              + "data: " + time_to_pretty_string(time) << endl;
  }

  if (collection.set_size != 0) {
    if (options.verbosity >= Verbosity::info) {
      cout << left << setw(report_col_width) << "Performance summary"
        << result.files.summary << endl;
      if (not options.evaluate.skip_viterbi_path)
        cout << left << setw(report_col_width) << "Viterbi path"
          << result.files.viterbi << endl;
      if (not options.evaluate.skip_bed)
        cout << left << setw(report_col_width) << "Motif occurrences (BED)"
          << result.files.bed << endl;
      if (not options.evaluate.skip_occurrence_table)
        cout << left << setw(report_col_width) << "Motif occurrences (table)"
          << result.files.table << endl;
    }

    if (not options.evaluate.skip_viterbi_path)
      viterbi_file.open(result.files.viterbi.c_str(), flags);
    if (not options.evaluate.skip_bed)
      bed_file.open(result.files.bed.c_str(), flags);
    if (not options.evaluate.skip_occurrence_table)
      occurrence_file.open(result.files.table.c_str(), flags);

    boost::iostreams::filtering_stream<boost::iostreams::output> v_out, bed_out,
        occ_out;
    switch (options.output_compression) {
      case Options::Compression::gzip:
        v_out.push(boost::iostreams::gzip_compressor());
        bed_out.push(boost::iostreams::gzip_compressor());
        occ_out.push(boost::iostreams::gzip_compressor());
        break;
      case Options::Compression::bzip2:
        v_out.push(boost::iostreams::bzip2_compressor());
        bed_out.push(boost::iostreams::bzip2_compressor());
        occ_out.push(boost::iostreams::bzip2_compressor());
        break;
      default:
        break;
    }
    v_out.push(viterbi_file);
    bed_out.push(bed_file);
    occ_out.push(occurrence_file);

    hmm.print_occurrence_table_header(occ_out);

    // TODO reactivate!
    if (false && not options.evaluate.skip_summary) {
      summary_out << endl;
      Training::Task my_task;
      my_task.measure = Measure::ClassificationPosterior;
      // TODO EVALUATION
      // TODO consider using compute_score instead of compute_score_all_motifs
      summary_out << "class log posterior = "
                  << hmm.compute_score_all_motifs(collection, my_task.measure,
                                                  options) << endl;
      my_task.measure = Measure::ClassificationLikelihood;
      // TODO EVALUATION
      // TODO consider using compute_score instead of compute_score_all_motifs
      summary_out << "class log likelihood = "
                  << hmm.compute_score_all_motifs(collection, my_task.measure,
                                                  options) << endl;
    }

    for (auto &contrast : collection) {
      summary_out << endl << endl << "Summary statistics for contrast '"
                  << contrast_name_tag(contrast) << "'" << endl;
      vector<ResultsCounts> counts;
      for (auto &dataset : contrast) {
        ResultsCounts c = evaluate_dataset(dataset, summary_out, v_out, occ_out,
                                           bed_out, options);
        counts.push_back(c);
      }

      if (not options.evaluate.skip_summary) {
        set<size_t> idxs;
        for (auto &a : counts)
          for (auto &b : a.exp_sites)
            idxs.insert(b.first);
        for (auto i : idxs) {
          if (hmm.is_motif_group(i)) {
            size_t motif_len = hmm.get_motif_len(i);
            string name = hmm.get_group_name(i);
            string consensus = hmm.get_group_consensus(i);

            string exp_motif_tag = "Expected motif counts - " + name
                                   + ":" + consensus;
            summary_out << endl
                        << "# Expected motif count = probabilistic count of "
                           "motif occurrences in the sequences" << endl
                        << exp_motif_tag << endl;
            matrix_t em(counts.size(), 2);
            for (size_t j = 0; j < counts.size(); j++) {
              em(j, 0) = counts[j].exp_motifs[i];
              em(j, 1) = contrast.sets[j].seq_size - em(j, 0);
            }
            count_report(summary_out, em, motif_len, contrast,
                         hmm.get_pseudo_count(), options.limit_logp, name,
                         tag + exp_motif_tag + " ");

            string vit_motif_tag = "Viterbi motif counts - " + name
                                   + ":" + consensus;
            summary_out << endl
                        << "# Viterbi motif count = number of motif "
                           "occurrences in the Viterbi paths of the sequences"
                        << endl
                        << vit_motif_tag << endl;
            matrix_t vm(counts.size(), 2);
            for (size_t j = 0; j < counts.size(); j++) {
              vm(j, 0) = counts[j].viterbi_motifs[i];
              vm(j, 1) = contrast.sets[j].seq_size - vm(j, 0);
            }
            count_report(summary_out, vm, motif_len, contrast,
                         hmm.get_pseudo_count(), options.limit_logp, name,
                         tag + vit_motif_tag + " ");

            string exp_tag = "Expected site counts - " + name + ":"
                             + consensus;
            summary_out << endl
                        << "# Expected site count = probabilistic count of "
                           "sequences with at least one motif occurrence" << endl
                        << exp_tag << endl;
            matrix_t e(counts.size(), 2);
            for (size_t j = 0; j < counts.size(); j++) {
              e(j, 0) = counts[j].exp_sites[i];
              e(j, 1) = contrast.sets[j].set_size - e(j, 0);
            }
            count_report(summary_out, e, motif_len, contrast,
                         hmm.get_pseudo_count(), options.limit_logp, name,
                         tag + exp_tag + " ");

            string vit_tag = "Viterbi site counts - " + name + ":"
                             + consensus;
            summary_out << endl
                        << "# Viterbi site count = number of sequences "
                           "for which the Viterbi path has at least one motif "
                           "occurrence" << endl
                        << vit_tag << endl;
            matrix_t v(counts.size(), 2);
            for (size_t j = 0; j < counts.size(); j++) {
              v(j, 0) = counts[j].viterbi_sites[i];
              v(j, 1) = contrast.sets[j].set_size - v(j, 0);
            }
            count_report(summary_out, v, motif_len, contrast,
                         hmm.get_pseudo_count(), options.limit_logp, name,
                         tag + vit_tag + " ");
          }
        }
      }
    }
  }
  return result;
}
