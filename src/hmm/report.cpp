
#include <iomanip>
#include <fstream>
#include <vector>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include "../aux.hpp"
#include "report.hpp"
#include "../timer.hpp"
#include "../plasma/find.hpp"
#include "../stats_config.hpp"

using namespace std;

void print_table(ostream &ofs, const matrix_t m, const Data::Series &data, size_t width, size_t prec) {
  size_t w = 0;
  for(size_t i = 0; i < m.size1(); i++) {
    size_t n = data.sets[i].path.size();
    if(n > w)
      w = n;
  }

  ofs << left << setw(w) << "" << right << setw(width) << "Present" << right << setw(width) << "Absent";
  if(m.size2() == 2)
    ofs << right << setw(width) << "Percent";
  ofs << endl;

  size_t prev_prec = ofs.precision();
  for(size_t i = 0; i < m.size1(); i++) {
    ofs << left << setw(w) << data.sets[i].path;
    for(size_t j = 0; j < m.size2(); j++)
      ofs << right << fixed << setprecision(prec) << setw(width) << m(i,j);
    if(m.size2() == 2)
      ofs << right << fixed << setprecision(prec) << setw(width) << (m(i,0) / (m(i,0) + m(i,1)) * 100);
    ofs << left << endl;
  }
  ofs.unsetf(ios_base::fixed);
  ofs.precision(prev_prec);
}

void count_report(ostream &ofs, const matrix_t counts, size_t motif_len, const Data::Series &data, double pseudo_count, bool limit_logp, const string &name, const string &prefix)
{
  size_t df = (counts.size1() - 1) * (counts.size2() - 1);
  print_table(ofs, counts, data, 10, 2);
  double mutualinf = calc_mutual_information(counts, pseudo_count, true, false, false);
  double exp_mutualinf = calc_mutual_information(counts, pseudo_count, true, true, false);
  double var_mutualinf = calc_mutual_information(counts, pseudo_count, true, true, true);
  double g = calc_g_test(counts, pseudo_count);
  // double correct_class = hmm.correct_classification(data);  // TODO: reactivate
  double log_p_g = pchisq(g, df, false, true);
  double cor_log_p_g_stringent = log(149) * motif_len + log_p_g;
  if(limit_logp)
    cor_log_p_g_stringent = min<double>(0, cor_log_p_g_stringent);

  ofs << prefix << "Discriminatory mutual information = " << mutualinf << " bit per sequence" << endl;
  ofs << prefix << "Expected discriminatory mutual information = " << exp_mutualinf << " bit per sequence" << endl;
  ofs << prefix << "Variance of discriminatory mutual information = " << var_mutualinf << " bit per sequence" << endl;
  ofs << prefix << "Std. dev. of discriminatory mutual information = " << sqrt(var_mutualinf) << " bit per sequence" << endl;
  ofs << prefix << "Z-score of discriminatory mutual information = " << exp_mutualinf / sqrt(var_mutualinf) << " bit per sequence" << endl;
  ofs << prefix << "G-test = " << g << endl;
  ofs << prefix << "Log-P(Chi-Square(G-Test)) = " << log_p_g << endl;
  ofs << prefix << "Bonferroni corrected log-P(Chi-Square(G-Test)) = " << cor_log_p_g_stringent << endl;
}

void eval_contrast(const HMM &hmm, const Data::Series &data, ostream &ofs, bool limit_logp, const string &tag)
{
  // double mi = hmm.mutual_information(data, contrast);
  // ofs << "Summed discriminatory mutual information = " << mi << " bit per sequence" << endl;

  for(size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
    if(hmm.is_motif_group(group_idx)) {
      size_t motif_len = hmm.get_motif_len(group_idx);
      vector_t v = hmm.posterior_atleast_one(data, group_idx);
      matrix_t counts(v.size(), 2);
      for(size_t i = 0; i < v.size(); i++) {
        counts(i,0) = v(i);
        counts(i,1) = data.sets[i].set_size - v(i);
      }
      double mcc = hmm.matthews_correlation_coefficient(data, group_idx); // TODO compute from count table
      // double llr = calc_log_likelihood_ratio(counts, hmm.pseudo_count);
      // double dips_tscore = hmm.dips_tscore(data, feature);
      double dips_sitescore = hmm.dips_sitescore(data, group_idx);
      // double correct_class = hmm.correct_classification(data);  // TODO: reactivate

      string name = hmm.get_group_name(group_idx);

      ofs << endl << tag << "Expected occurrence statistics for motif " << name << endl;
      count_report(ofs, counts, motif_len, data, hmm.get_pseudo_count(), limit_logp, name, tag);
      ofs << tag << "Matthews correlation coefficient = " << mcc << endl;
      // ofs << "DIPS t-score = " << dips_tscore * 100 << " %" << endl;
      ofs << tag << "DIPS site-score = " << dips_sitescore * 100 << " %" << endl;
      // ofs << "log P correct classification = " << correct_class << endl; // TODO: reactivate
      // TODO print a table of the likelihoods
      ofs << tag << "Log likelihood difference = " << hmm.log_likelihood_difference(data, group_idx) << endl; // TODO compute from count table
    }
}

template <class X, class Y> double cor_pearson(const std::vector<X> &x, const std::vector<Y> &y)
{
  size_t n = x.size();
  double zx = std::accumulate(x.begin(), x.end(), 0);
  double zy = std::accumulate(y.begin(), y.end(), 0);
  zx /= n;
  zy /= n;

  double num = 0;
  double denom_x = 0, denom_y = 0;
  for(size_t i = 0; i < n; i++) {
    double a = x[i] - zx;
    double b = y[i] - zy;
    num += a*b;
    denom_x += a*a;
    denom_y += b*b;
  }
  if(num == 0)
    return(0);
  double r = num / sqrt(denom_x) / sqrt(denom_y);
  return(r);
}

template <class X> std::vector<double> tied_ranking(const std::vector<X> &x)
{
  const size_t n = x.size();
  std::vector<double> r(n);
//  std::iota(r.begin(), r.end(), 0);
  std::map<X,size_t> counts;
  for(auto &y: x)
    counts[y]++;
  size_t cumul = 0;
  for(auto &y: counts) {
    double first = cumul;
    cumul += y.second;
    double ave = (first+cumul-1) / 2;
    y.second = ave;
  }
  for(size_t i = 0; i < n; i++)
    r[i] = counts[x[i]];
  return(r);
}

template <class X, class Y> double cor_spearman(const std::vector<X> &x, const std::vector<Y> &y)
{
  return(cor_pearson(tied_ranking(x),tied_ranking(y)));
}

template <class X> double cor_rank(const std::vector<X> &x)
{
  std::vector<X> y(x.size());
  std::iota(y.begin(), y.end(), 0);
  return(cor_pearson(tied_ranking(x),y));
}

double cor_fisher_z(double x, size_t n)
{
  double z = sqrt((n-3) / 1.06) * atanh(x);
  return(z);
}

double cor_student_t(double x, size_t n)
{
  double t = x * sqrt((n-2) / ( 1-x*x));
  return(t);
}

string gen_stars(double x) {
  if(x < 0.001)
    return("***");
  if(x < 0.01)
    return("**");
  if(x < 0.05)
    return("*");
  return("");
}

const boost::math::normal_distribution<double> standard_normal_distribution;

template <class T> void correlation_report(const std::vector<T> &x, std::ostream &out, size_t width=12, size_t prec=5)
{
  double rho = cor_rank(x);
  size_t n = x.size();
  double z = cor_fisher_z(rho, n);
  double t = cor_student_t(rho, n);
  // out << "Rank correlation: rho = " << rho << " z = " << z << " t = " << t << std::endl;
  // using namespace boost::math;
  const boost::math::students_t_distribution<> students_t_dist(n-2);
  double p_norm = 0, p_t = 0;
  if(rho < 0) {
    p_t = boost::math::cdf(students_t_dist, t);
    p_norm = boost::math::cdf(standard_normal_distribution, z);
  } else {
    p_t = boost::math::cdf(complement(students_t_dist, t));
    p_norm = boost::math::cdf(complement(standard_normal_distribution, z));
  }
  size_t prev_prec = out.precision();
  out << setw(width) << fixed << right << setprecision(prec) << rho
    << setw(width) << fixed << right << setprecision(2) << z
    << setw(width) << fixed << right << setprecision(2) << log(p_norm)
    << setw(4) << right << gen_stars(p_norm)
    << setw(width) << fixed << right << setprecision(2) << t
    << setw(width) << fixed << right << setprecision(2) << log(p_t)
    << setw(4) << right << gen_stars(p_t)
    << std::endl;
  out.unsetf(ios_base::fixed);
  out.precision(prev_prec);
}

ResultsCounts evaluate_hmm_single_data_set(const HMM &hmm,
    const Data::Set &data,
    ostream &out,
    ostream &v_out,
    ostream &occurrence_out,
    const hmm_options &options)
{
  const size_t width = 12;
  const size_t prec = 5;
  Timer timer;

  map<size_t,double> n_sites;
  map<size_t,double> n_motifs;
  map<size_t,size_t> n_viterbi_sites;
  map<size_t,size_t> n_viterbi_motifs;

  if(options.evaluate.summary) {
    double log_likelihood = hmm.log_likelihood(data);
    // out << endl << "Summary of " << data.path << endl;
    // out << "Total sequences = " << data.sequences.size() << endl;
    out << "Log-likelihood of " << data.path << " = " << log_likelihood << endl;
    // out << "Akaike information criterion AIC = " << 2 * hmm.n_parameters() - 2 * log_likelihood << endl;
    //      out << "Akaike information criterion AIC = " << 2 * hmm.non_zero_parameters(hmm.gen_training_targets(options)) - 2 * log_likelihood << endl;
  }

  if(options.evaluate.viterbi_path)
    v_out << "# " << data.path << " details following" << endl;

  const size_t n_groups = hmm.get_ngroups();
  const size_t n = data.sequences.size();
  size_t number_motifs = 0;
  for(size_t group_idx = 0; group_idx < n_groups; group_idx++)
    if(hmm.is_motif_group(group_idx))
      number_motifs++;

  vector<vector<double>> atl_counts(number_motifs), exp_counts(number_motifs), vit_counts(number_motifs);
  for(size_t group_idx = 0; group_idx < number_motifs; group_idx++) {
    atl_counts[group_idx] = vector<double>(n);
    exp_counts[group_idx] = vector<double>(n);
    vit_counts[group_idx] = vector<double>(n);
  }

  if(options.evaluate.ric)
    for(size_t group_idx = 0; group_idx < n_groups; group_idx++)
      if(hmm.is_motif_group(group_idx)) {
        double ric = hmm.rank_information(data, group_idx);
        out << "RIC = " << ric << std::endl;
      }

  for(size_t i = 0; i < data.sequences.size(); i++) {
    HMM::StatePath path;
    double lp = hmm.viterbi(data.sequences[i], path);

    if(options.evaluate.viterbi_path or options.evaluate.summary) {
      stringstream viterbi_str, exp_str, atl_str;

      bool first = true;
      size_t motif_idx = 0;
      for(size_t group_idx = 0; group_idx < n_groups; group_idx++)
        if(hmm.is_motif_group(group_idx)) {
          if(first) first = false; else { viterbi_str << "/"; exp_str << "/"; atl_str << "/"; }

          double atl = hmm.posterior_atleast_one(data.sequences[i], group_idx).posterior;
          double expected = hmm.expected_posterior(data.sequences[i], group_idx);
          size_t n_viterbi = hmm.count_motif(path, group_idx);
          atl_counts[motif_idx][i] = atl;
          exp_counts[motif_idx][i] = expected;
          vit_counts[motif_idx][i] = n_viterbi;

          n_sites[group_idx] += atl;
          n_motifs[group_idx] += expected;
          n_viterbi_sites[group_idx] += (n_viterbi > 0 ? 1 : 0);
          n_viterbi_motifs[group_idx] += n_viterbi;

          if(options.evaluate.viterbi_path) {
            viterbi_str << n_viterbi;
            exp_str << expected;
            atl_str << atl;
          }

          motif_idx++;
        }
      if(options.evaluate.viterbi_path) {
        v_out << ">" << data.sequences[i].definition << endl;
        v_out << "V-sites = " << viterbi_str.str() << " E-sites = " << exp_str.str() << " P(#sites>=1) = " << atl_str.str() << " Viterbi log-p = " << lp << endl;
        v_out << data.sequences[i].sequence << endl;
        v_out << hmm.path2string_group(path) << endl;
      }
    }

//    if(options.print_posterior) {
//      vector_t posterior = hmm.motif_posterior(data.sequences[i]);
//      v_out << "Posterior ";
//      for(auto v : posterior)
//        v_out << " " << v;
//      v_out << endl;
//    }

    if(options.evaluate.occurrence_table) {
      // print to motif occurrence table
      for(size_t pos = 0; pos < path.size(); pos++)
        for(size_t group_idx = 0; group_idx < n_groups; group_idx++)
          if(hmm.is_motif_group(group_idx) and
              path[pos] == hmm.groups[group_idx].states[0]) {
            size_t end = pos + 1;
            while(end != path.size() and path[end] > path[pos])
              end++;
            string seq = data.sequences[i].sequence.substr(pos, end-pos);
            occurrence_out << data.path << "\t" << data.sequences[i].definition << "\t" << pos << "\t" << hmm.groups[group_idx].name << "\t" << seq << endl;
          }
    }
  }

  if(options.evaluate.summary) {
    out << setw(30) << left << "Rank analysis";
    out << setw(width) << right << "Rho"
      << setw(width) << right << "Z"
      << setw(width) << right << "log P(Z)"
      << setw(4) << right << ""
      << setw(width) << right << "t"
      << setw(width) << right << "log P(t)"
      << setw(4) << right << "" << std::endl;
    // TODO: write out the motif name
    for(size_t group_idx = 0; group_idx < number_motifs; group_idx++) {
      out << setw(30) << left << "Posterior decoded site";
      correlation_report(atl_counts[group_idx], out, width, prec);
      out << setw(30) << left << "Posterior number motif";
      correlation_report(exp_counts[group_idx], out, width, prec);
      vector<size_t> vit(n);
      for(size_t i = 0; i < n; i++)
        vit[i] = vit_counts[group_idx][i] > 0;
      out << setw(30) << left << "Viterbi site";
      correlation_report(vit, out, width, prec);
      out << setw(30) << left << "Viterbi number motif";
      correlation_report(vit_counts[group_idx], out, width, prec);
    }
  }

  double time = timer.tock();
  if(options.timing_information)
    cerr << "Evaluation for " << data.path << ": " << time << " micro-seconds" << endl;
  ResultsCounts results = {n_sites, n_motifs, n_viterbi_sites, n_viterbi_motifs};
  return(results);
}

void evaluate_hmm(const HMM &hmm,
    const Data::Collection &data,
    const string &tag,
    const Training::Tasks &tasks,
    const hmm_options &options)
{
  /*  TODO: re-enable class-based models
  if(options.adapt_classes) {
    if(options.verbosity >= Verbosity::info)
      cout << "Registering data sets for class based HMMs." << endl;

    for(auto &series: data)
      for(auto &data: series)
        hmm.register_class_hmm(data);
    }
    Training::Targets generative_targets = hmm.gen_secondary_training_targets(options.training_type);
    Training::Task generative_task;
    generative_task.method = Training::Method::reestimation;
    generative_task.measure = Measure::likelihood;
    generative_task.targets = generative_targets;
    hmm.reestimation(data, generative_task, options);
  } */

  ios_base::openmode flags = ios_base::out;
  if(options.output_compression != Compression::none)
    flags |= ios_base::binary;

  string file_tag = "";
  if(tag != "")
    file_tag = tag + ".";
  string eval_out_path = options.label + file_tag + ".summary";
  string occurrence_out_path = options.label + file_tag + ".table" + compression2ending(options.output_compression);
  string viterbi_out_path = options.label + file_tag + ".viterbi" + compression2ending(options.output_compression);

  ofstream summary_out, occurrence_file, viterbi_file;
  summary_out.open(eval_out_path.c_str());

  if(options.evaluate.summary) {
    size_t col0_w = 10, col1_w = 10, col2_w = 10;
    for(size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
      if(hmm.is_motif_group(group_idx)) { // do not output information about non-motif groups
        string consensus = hmm.get_group_consensus(group_idx);
        string name = hmm.get_group_name(group_idx);
        if(consensus.size() > col1_w)
          col1_w = consensus.size();
        if(name.size() > col0_w)
          col0_w = name.size();
      }
    col0_w += 2; col1_w += 2;
    summary_out << "Motif summary" << endl;
    // TODO: print out PWM and affinity matrix
    summary_out << left << setw(col0_w) << "Motif name"
      << left << setw(col1_w) << "Consensus"
      << right << setw(col2_w) << "IC [bit]" << endl;
    for(size_t group_idx = 0; group_idx < hmm.get_ngroups(); group_idx++)
      if(hmm.is_motif_group(group_idx)) { // do not output information about non-motif groups
        string consensus = hmm.get_group_consensus(group_idx);
        string name = hmm.get_group_name(group_idx);
        double ic = hmm.information_content(group_idx);
        summary_out << left << setw(col0_w) << name
          << left << setw(col1_w) << consensus
          << right << setw(col2_w) << ic << endl;
      }
    summary_out << endl;

    Timer eval_timer;
    for(auto &series: data)
      eval_contrast(hmm, series, summary_out, options.limit_logp, tag);
    double time = eval_timer.tock();

    if(options.timing_information)
      cerr << "Evaluation of contrast for " << (tag == "" ? "" : tag + " ") << "data: " << time << " micro-seconds" << endl;
  }

  if(data.set_size != 0) {
    if(options.verbosity >= Verbosity::info)
      cout << "Performance summary in " << eval_out_path << endl
        << "Viterbi path in " << viterbi_out_path << endl
        << "Motif occurrence table in " << occurrence_out_path << endl;

    occurrence_file.open(occurrence_out_path.c_str());
    viterbi_file.open(viterbi_out_path.c_str(), flags);

    boost::iostreams::filtering_stream<boost::iostreams::output> v_out, occurrence_out;
    switch(options.output_compression) {
      case Compression::gzip:
        v_out.push(boost::iostreams::gzip_compressor());
        occurrence_out.push(boost::iostreams::gzip_compressor());
        break;
      case Compression::bzip2:
        v_out.push(boost::iostreams::bzip2_compressor());
        occurrence_out.push(boost::iostreams::bzip2_compressor());
        break;
      default:
        break;
    }
    v_out.push(viterbi_file);
    occurrence_out.push(occurrence_file);


    if(options.evaluate.summary) {
      summary_out << endl;
      Training::Task my_task;
      my_task.measure = Measure::ClassificationPosterior;
      // TODO consider using compute_score instead of compute_score_all_motifs
      summary_out << "class log posterior = " << hmm.compute_score_all_motifs(data, my_task.measure, options.weighting) << std::endl;
      my_task.measure = Measure::ClassificationLikelihood;
      // TODO consider using compute_score instead of compute_score_all_motifs
      summary_out << "class log likelihood = " << hmm.compute_score_all_motifs(data, my_task.measure, options.weighting) << std::endl;
    }

    for(auto &series: data) {
      vector<ResultsCounts> counts;
      for(auto &data_set: series) {
        ResultsCounts c = evaluate_hmm_single_data_set(hmm, data_set, summary_out, v_out, occurrence_out, options);
        counts.push_back(c);
      }

      if(options.evaluate.summary) {
        set<size_t> idxs;
        for(auto &a: counts)
          for(auto &b: a.exp_sites)
            idxs.insert(b.first);
        for(auto i: idxs) {
          if(hmm.is_motif_group(i)) {
            size_t motif_len = hmm.get_motif_len(i);
            string name = hmm.get_group_name(i);

            string exp_motif_tag = "Posterior decoded motif counts - " + name;
            summary_out << endl << exp_motif_tag << endl;
            matrix_t em(counts.size(),2);
            for(size_t j = 0; j < counts.size(); j++) {
              em(j,0) = counts[j].exp_motifs[i];
              em(j,1) = series.sets[j].seq_size - em(j,0);
            }
            count_report(summary_out, em, motif_len, series, hmm.get_pseudo_count(), options.limit_logp, name, tag + exp_motif_tag + " ");

            string vit_motif_tag = "Viterbi decoded motif counts - " + name;
            summary_out << endl << vit_motif_tag << endl;
            matrix_t vm(counts.size(),2);
            for(size_t j = 0; j < counts.size(); j++) {
              vm(j,0) = counts[j].viterbi_motifs[i];
              vm(j,1) = series.sets[j].seq_size - vm(j,0);
            }
            count_report(summary_out, vm, motif_len, series, hmm.get_pseudo_count(), options.limit_logp, name, tag + vit_motif_tag + " ");


            string exp_tag = "Posterior decoded site counts - " + name;
            summary_out << endl << exp_tag << endl;
            matrix_t e(counts.size(),2);
            for(size_t j = 0; j < counts.size(); j++) {
              e(j,0) = counts[j].exp_sites[i];
              e(j,1) = series.sets[j].set_size - e(j,0);
            }
            count_report(summary_out, e, motif_len, series, hmm.get_pseudo_count(), options.limit_logp, name, tag + exp_tag + " ");

            string vit_tag = "Viterbi decoded site counts - " + name;
            summary_out << endl << vit_tag << endl;
            matrix_t v(counts.size(),2);
            for(size_t j = 0; j < counts.size(); j++) {
              v(j,0) = counts[j].viterbi_sites[i];
              v(j,1) = series.sets[j].set_size - v(j,0);
            }
            count_report(summary_out, v, motif_len, series, hmm.get_pseudo_count(), options.limit_logp, name, tag + vit_tag + " ");
          }
        }
      }
    }
  }
}
 
