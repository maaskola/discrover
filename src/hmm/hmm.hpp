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
 *       Filename:  hmm.hpp
 *
 *    Description:  Code for binding site hidden Markov models
 *
 *        Created:  Wed Aug 3 02:08:55 2011 +0200
 *
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef HMM_HPP
#define HMM_HPP

/* Description

   The code in this package allows to train Hidden Markov Models using different learning paradigms.

   1. Signal only learning
     1.a. Maximum likelihood learning
       - Baum-Welch, scaled and unscaled
       - Likelihood gradient, scaled
     1.b. Viterbi learning
     1.c. Posterior learning
       - Gradient baseed, scaled; This does not actually work in practice but is required for 2.a.
   2. Discriminative learning
     2.a. Discriminatory mutual information
       - Gradient learning, scaled
     2.b. Difference of likelihood in signal and control, like in the program DME
       - Gradient learning, scaled
     2.c. Difference of motif occurrence probability between signal and control, like in the program DIPS -- TODO
       - Gradient learning, scaled
     2.d. Difference of site occurrence probability between signal and control, similar to the program DIPS
       - Gradient learning, scaled
     2.e. Log of the joint probability of classifying all samples correctly
       - Gradient learning, scaled
     2.f. ML learning of the full model on the signal, and a reduced model on the control data
       - Modified Baum-Welch, scaled
     2.g. Matthew's correlation coefficient, a.k.a. Phi coefficient, a.k.a. Cramer's V
       - Gradient learning, scaled

   With all gradient based learning approaches it is possible to learn both transition and emission
   probabilities, or just one of the two sets of probabilities.
*/


// TODO
// Discriminative learning:
// * D  implement topological extension for conditional modeling
// * E2  implement separation based on expected posterior motif probability gradient
//
// Test these approaches' generalization capabilities:
//   ML-based, signal only:
//      Baum-Welch on signal only
//      Viterbi on signal only
//      Likelihood gradient on signal only
//   ML-based, discriminative
//      Baum-Welch on extended topology
//      Likelihood gradient on extended topology
//   MMI-based
//      MI gradient
//

// #define BOOST_UBLAS_MOVE_SEMANTICS // TODO find out if this helps

#include <boost/container/map.hpp>
#include <boost/container/flat_map.hpp>
#include <list>
#include <unordered_map>
#include "association.hpp"
#include "hmm_options.hpp"
#include "../verbosity.hpp"


struct Gradient {
  Gradient() : transition(), emission() { };
  matrix_t transition;
  matrix_t emission;
};

inline double scalar_product(const Gradient &gradient1, const Gradient &gradient2)
{
  assert(gradient1.transition.size1() == gradient2.transition.size1());
  assert(gradient1.transition.size2() == gradient2.transition.size2());
  double g = 0;
  for(size_t i = 0; i < gradient1.transition.size1(); i++)
    for(size_t j = 0; j < gradient1.transition.size2(); j++)
      g += gradient1.transition(i,j) * gradient2.transition(i,j);
  for(size_t i = 0; i < gradient1.emission.size1(); i++)
    for(size_t j = 0; j < gradient1.emission.size2(); j++)
      g += gradient1.emission(i,j) * gradient2.emission(i,j);
  return(g);
}

/** The directional derivate */
inline double dderiv(const Gradient &direction, const Gradient &gradient)
{
  assert(fabs(sqrt(scalar_product(direction,direction)) - 1) < 1e-6);
  double x = scalar_product(direction, gradient); // / sqrt(scalar_product(direction, direction));
  return(x);
}


class HMM;

struct ResultsCounts;

ResultsCounts evaluate_hmm_single_data_set(const HMM &hmm, const Data::Set &data, std::ostream &out, std::ostream &v_out, std::ostream &occurrence_out, const hmm_options &options);
double train_hmm(HMM &hmm, const Data::Collection &training_data, const Training::Tasks &tasks, const hmm_options &options);

struct Group {
  enum class Kind {
    Special,
    Background,
    Motif
  };
  Kind kind;
  std::string name;
  Training::Range states;
};

class HMM {
  public:
    // constructors
    HMM(size_t bg_order, Verbosity verbosity_, double pseudo_count_=1.0);
    HMM(const std::string &path, Verbosity verbosity_, double pseudo_count_=1.0);
    HMM(const HMM &hmm, bool copy_deep=true);

    typedef std::unordered_map<std::string, std::vector<size_t>> mask_sub_t;
    typedef std::unordered_map<std::string, mask_sub_t> mask_t;
    typedef boost::numeric::ublas::vector<size_t> StatePath;

  protected:
    /** The size of the alphabet. */
    const static size_t alphabet_size = 4;
    /** The number of emissions. */
    size_t n_emissions; // = sum(alphabet_size**i, i, 1, bg_order+1)
    /** Level of output verbosity. */
    Verbosity verbosity;
    /** Whether to save intermediate parameters on disc during learning. */
    bool store_intermediate;
    /** The index of the start state. */
    const static size_t start_state = 0;
    /** The index of the bg state. */
    const static size_t bg_state = 1;
    /** The index of the first motif state. */
    const static size_t first_state = 2;
    /** The index of the last motif state. */
    size_t last_state;
    /** The number of states. */
    size_t n_states;
    /** The pseudo count to add to contingency tables. */
    double pseudo_count;

    std::vector<Group> groups;
    /** Classification of each node */
    std::vector<size_t> group_ids;
    std::vector<size_t> order;
    std::vector<size_t> order_offset;
    size_t max_order;
    const static size_t initial_history = 0;

    /** The transition probabilities. */
    matrix_t transition;
    /** The emission probabilities. */
    matrix_t emission;

    /** The indices of the predecessors of each state. */
    std::vector<std::list<size_t>> pred;
    /** The indices of the successors of each state. */
    std::vector<std::list<size_t>> succ;

    Training::Range all_range, constitutive_range, emitting_range;


    struct RegisteredDataSet {
      Specification::DataSet spec;
      double class_prior;
      std::map<size_t, double> motif_prior; // the key is supposed to be the group index of the motif
      double get_motif_prior(size_t group_idx) const;
    };
    std::unordered_map<std::string, RegisteredDataSet> registered_datasets; // the keys are SHA1 hashes of the file contents
    double compute_marginal_motif_prior(size_t group_idx) const;
    double get_class_motif_prior(const std::string &sha1, size_t group_idx) const;
    double get_class_prior(const std::string &sha1) const;

// -------------------------------------------------------------------------------------------
// Initialization routines
// -------------------------------------------------------------------------------------------

    void initialize_emissions(size_t bg_order);
    void initialize_transitions();
    void initialize_bg_transitions();
    void initialize_order_offsets();
    void initialize_transitions_to_and_from_chain(size_t w, double l, double lambda, size_t first, size_t last, size_t pad_left, size_t pad_right);
    void set_motif_emissions(const matrix_t &e, size_t first, size_t n_insertions, size_t pad_left, size_t pad_right);
    void normalize_transition(matrix_t &m) const;
    void normalize_emission(matrix_t &m) const;

    /** Initialize the predecessor and successor data structures */
    void initialize_pred_succ();

    /** Initialize the range data structures */
    void initialize_ranges();

    /** Initialize the predecessor&successor and range data structures */
    void finalize_initialization();

    /**  Check whether the motif is enriched in the desired samples */
    bool check_enrichment(const Data::Series &data, const matrix_t &counts, size_t group_idx) const;

  public:
    /* Add a motif based on a IUPAC string. */
    size_t add_motif(const std::string &seq, double alpha, double exp_seq_len, double lambda, const std::string &name, const std::vector<size_t> &insertions, size_t pad_left, size_t pad_right);
    /* Add a motif based on a matrix. */
    size_t add_motif(const matrix_t &e, double exp_seq_len, double lambda, const std::string &name, std::vector<size_t> insertions, size_t pad_left, size_t pad_right);

    /* Add motifs of another HMM. */
    void add_motifs(const HMM &hmm);

    /** Perform HMM training. */
    double train(const Data::Collection &data, const Training::Tasks &tasks, const hmm_options &options);

    /** Perform iterative HMM training, like re-estimation (expectation maximization) or gradient based. */
    void iterative_training(const Data::Collection &data, const Training::Tasks &tasks, const hmm_options &options);
    /** Perform one iteration of iterative HMM training. */
    bool perform_training_iteration(const Data::Collection &data, const Training::Tasks &tasks, const hmm_options &options, Training::State &ts);
    /** Perform one iteration of gradient training. */
    bool perform_training_iteration_gradient(const Data::Collection &data, const Training::Task &task, const hmm_options &options, int &center, double &score);
    /** Perform one iteration of re-estimation training. */
    bool perform_training_iteration_reestimation(const Data::Collection &data, const Training::Task &task, const hmm_options &options, double &score);
    /** Perform one iteration of class parameter re-estimation training. */
    bool reestimate_class_parameters(const Data::Collection &data, const Training::Task &task, const hmm_options &options, double &score);

    /** Perform HMM training for the background part of the model using Baum-Welch. */
    void train_background(const Data::Collection &data, const hmm_options &options);

// -------------------------------------------------------------------------------------------
// Saving and loading
// -------------------------------------------------------------------------------------------

    /** Format parameters in a parseable text-format. */
    void serialize(std::ostream &os, const ExecutionInformation &exec_info, size_t format_version=6) const;
    /** Restore parameters from parseable text-format. */
    void deserialize(std::istream &os);

// -------------------------------------------------------------------------------------------
// Some auxiliary routines
// -------------------------------------------------------------------------------------------

    /** The emission information content of the HMM. */
    double information_content(size_t motif_idx) const;

    /** A routine that calculates the emission consensus. */
    std::string get_group_consensus(size_t idx, double threshold=0.1) const;

    /** A getter routine for the group names. */
    std::string get_group_name(size_t idx) const;

    /** A getter routine for the number of groups, including the constitutive one. */
    size_t get_ngroups() const;

    /** A getter routine for the number of motifs, excluding the constitutive group. */
    size_t get_nmotifs() const;

    /** A getter routine for the length of a motif. */
    size_t get_motif_len(size_t motif_idx) const;

    /** A getter routine for the pseudo count. */
    double get_pseudo_count() const;

    /** Generate a random candidate variant for MCMC sampling */
    HMM random_variant(const hmm_options &options, std::mt19937 &rng) const;

    size_t count_motif(const StatePath &path, size_t motif) const;

    /** Prepare training tasks according to the objectives given in the options. */
    Training::Tasks define_training_tasks(const hmm_options &options) const;

    size_t n_parameters() const;
    size_t non_zero_parameters(const Training::Targets &targets) const;

    friend std::ostream &operator<<(std::ostream& os, const HMM &hmm);
    friend ResultsCounts evaluate_hmm_single_data_set(const HMM &hmm, const Data::Set &data, std::ostream &out, std::ostream &v_out, std::ostream &occurrence_out, const hmm_options &options);
    friend double train_hmm(HMM &hmm, const Data::Collection &training_data, const Training::Tasks &tasks, const hmm_options &options);

    double compute_score(const Data::Collection &data, const Training::Task &task) const;

    void shift_forward(size_t group_idx, size_t n);
    void shift_backward(size_t group_idx, size_t n);

// -------------------------------------------------------------------------------------------
// Probabilistic evaluation
// -------------------------------------------------------------------------------------------

    struct posterior_t {
      double log_likelihood;
      double posterior;
    };
    struct posterior_gradient_t {
      double log_likelihood;
      double posterior;
      matrix_t T;
      matrix_t E;
    };

    vector_t posterior_atleast_one(const Data::Series &data, size_t group_idx) const;
    vector_t posterior_atleast_one(const Data::Set &data, size_t group_idx) const;
    double   sum_posterior_atleast_one(const Data::Set &data, size_t group_idx) const;
    posterior_t posterior_atleast_one(const Data::Seq &data, size_t group_idx) const;

    vector_t viterbi_atleast_one(const Data::Series &data, size_t group_idx) const;
    double   viterbi_atleast_one(const Data::Set &data, size_t group_idx) const;

    double viterbi(const Data::Seq &s, StatePath &path) const;
    double viterbi_zeroth_order(const Data::Seq &s, StatePath &path) const;
    double viterbi_higher_order(const Data::Seq &s, StatePath &path) const;

    vector_t expected_posterior(const Data::Series &data, size_t group_idx) const;
    double   expected_posterior(const Data::Set &data, size_t group_idx) const;
    double   expected_posterior(const Data::Seq &data, size_t group_idx) const;

// -------------------------------------------------------------------------------------------
// Generative and discriminative measures
// -------------------------------------------------------------------------------------------

    double log_likelihood(const Data::Collection &data) const;
    double log_likelihood(const Data::Series &data) const;
    double log_likelihood(const Data::Set &s) const;

    double class_likelihood(const Data::Series &data, const std::vector<size_t> &groups, bool compute_posterior) const;
    double class_likelihood(const Data::Series &data, size_t group_idx, bool compute_posterior) const;
    double class_likelihood(const Data::Set &data, size_t group_idx, bool compute_posterior) const;

    double log_likelihood_difference(const Data::Series &data, const std::vector<size_t> &groups) const;
    double log_likelihood_difference(const Data::Series &data, size_t group_idx) const;
    /** The mutual information of condition and motif occurrence */
    double mutual_information(const Data::Series &data, const std::vector<size_t> &groups) const;
    /** The mutual information of condition and motif occurrence */
    double mutual_information(const Data::Series &data, size_t group_idx) const;
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Series &data, const std::vector<size_t> &groups) const;
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Series &data, size_t group_idx) const;
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Set &data, size_t group_idx) const;
    /** The summed Matthew's correlation coefficients of all individual contrasts. */
    double matthews_correlation_coefficient(const Data::Series &data, const std::vector<size_t> &groups) const;
    /** The summed Matthew's correlation coefficients of all individual contrasts. */
    double matthews_correlation_coefficient(const Data::Series &data, size_t group_idx) const;
    /** The summed difference of expected occurrence frequencies. */
    double dips_tscore(const Data::Series &data, const std::vector<size_t> &groups) const;
    /** The summed difference of expected occurrence frequencies. */
    double dips_tscore(const Data::Series &data, size_t group_idx) const;
    /** The summed difference of site occurrence. */
    double dips_sitescore(const Data::Series &data, const std::vector<size_t> &groups) const;
    /** The summed difference of site occurrence. */
    double dips_sitescore(const Data::Series &data, size_t group_idx) const;

// -------------------------------------------------------------------------------------------
// Expectation-maximization type learning
// -------------------------------------------------------------------------------------------

    /** Perform full expectation-maximization type learning */
    void reestimation(const Data::Collection &data, const Training::Task &task, const hmm_options &options);

  protected:
    /** Perform one iteration of expectation-maximization type learning */
    double reestimationIteration(const Data::Collection &data, const Training::Task &task, const hmm_options &options);

    void M_Step(const matrix_t &T, const matrix_t &E, const Training::Targets &targets, const hmm_options &options);

    /** Perform one iteration of scaled Baum-Welch learning for a data collection */
    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Collection &data, const Training::Targets &targets, const hmm_options &options) const;
    /** Perform one iteration of scaled Baum-Welch learning for a data set */
    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Set &data, const Training::Targets &targets, const hmm_options &options) const;
    /** Perform one iteration of Viterbi learning for a data collection */
    double ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Collection &data, const Training::Targets &targets, const hmm_options &options);
    /** Perform one iteration of Viterbi learning for a data set */
    double ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Set &data, const Training::Targets &targets, const hmm_options &options);

    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const;
    double BaumWelchIteration_single(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const;

// -------------------------------------------------------------------------------------------
// Monte-Carlo Markov Chain inference
// -------------------------------------------------------------------------------------------
//
    /** Perform parallel tempering */
    std::vector<std::list<std::pair<HMM, double>>> mcmc(const Data::Collection &data, const Training::Task &task, const hmm_options &options);

// -------------------------------------------------------------------------------------------
// Gradient learning
// -------------------------------------------------------------------------------------------

    /** Perform gradient line searching with the desired objective function
     * using an approximate exponential step search.  **/
    std::pair<double, HMM> line_search(const Data::Collection &data, const Gradient &gradient, const Training::Task &task, int &center) const;
    /** Perform gradient line searching with the Moré-Thuente algorithm and the
     * desired objective function. */
    std::pair<double, HMM> line_search_more_thuente(const Data::Collection &data, const Gradient &gradient, double score, int &info, const Training::Task &task, const hmm_options &options) const;
    /** Auxiliary routine to build candidate HMM for a given direction and step
     * size. */
    HMM build_trial_model(const Gradient &gradient, double alpha, const Training::Task &task) const;

    /** Compute the gradient of the desired objective function */
    Gradient compute_gradient(const Data::Collection &data, double &score, const Training::Task &task) const;
    Gradient compute_gradient(const Data::Series &data, double &score, const Training::Task &task) const;
    double compute_gradient(const Data::Series &data, Gradient &gradient, const Training::Task &task, size_t group_idx) const;


// -------------------------------------------------------------------------------------------
// Forward and backward algorithms, and associated likelihood and posterior calculation
// -------------------------------------------------------------------------------------------

    /** The standard forward algorithm with scaling.
     *  The scaling vector is also determined. */
    matrix_t compute_forward_scaled(const Data::Seq &s, vector_t &scale) const;
    matrix_t compute_forward_scaled_zeroth_order(const Data::Seq &s, vector_t &scale) const;
    matrix_t compute_forward_scaled_higher_order(const Data::Seq &s, vector_t &scale) const;
    /** Computes only the scaling vector of the standard forward algorithm with scaling. */
    vector_t compute_forward_scale(const Data::Seq &s) const;
    vector_t compute_forward_scale_zeroth_order(const Data::Seq &s) const;
    vector_t compute_forward_scale_higher_order(const Data::Seq &s) const;

    /** The standard forward algorithm with prescaling.
     *  The scaling vector is assumed to be given. */
    matrix_t compute_forward_prescaled(const Data::Seq &s, const vector_t &scale) const;
    matrix_t compute_forward_prescaled_zeroth_order(const Data::Seq &s, const vector_t &scale) const;
    matrix_t compute_forward_prescaled_higher_order(const Data::Seq &s, const vector_t &scale) const;
    /** The standard backward algorithm with prescaling.
     *  The scaling vector is assumed to be given. */
    matrix_t compute_backward_prescaled(const Data::Seq &s, const vector_t &scale) const;
    matrix_t compute_backward_prescaled_zeroth_order(const Data::Seq &s, const vector_t &scale) const;
    matrix_t compute_backward_prescaled_higher_order(const Data::Seq &s, const vector_t &scale) const;

    double likelihood_from_scale(const vector_t &scale) const;
    double log_likelihood_from_scale(const vector_t &scale) const;

    double expected_state_posterior(size_t k, const matrix_t &f, const matrix_t &b, const vector_t &scale) const;


// -------------------------------------------------------------------------------------------
// Gradient methods
//
// These methods compute for various functions the gradient with respect to transformed
// transition and emission probabilities of the HMM.
// -------------------------------------------------------------------------------------------

    /** Likelihood gradient with respect to transformed transition and emission probabilities */
    double log_likelihood_gradient(const Data::Seqs &s, const Training::Targets &targets, matrix_t &transition_g, matrix_t &emission_g) const;

    /** Gradient of mutual information of class and motif occurrence with respect to transformed transition and emission probabilities */
    double mutual_information_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of mutual information of rank and motif occurrence with respect to transformed transition and emission probabilities */
    double rank_information_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;
    double rank_information_gradient(const Data::Set &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of Matthew's correlation coefficient of class and motif occurrence with respect to transformed transition and emission probabilities */
    double matthews_correlation_coefficient_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of Chi-Square of class and motif occurrence with respect to transformed transition and emission probabilities */
    double chi_square_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of log likelihood difference with respect to transformed transition and emission probabilities */
    double log_likelihood_difference_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of site frequency difference with respect to transformed transition and emission probabilities */
    double site_frequency_difference_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of class likelihood with respect to transformed transition and emission probabilities */
    double class_likelihood_gradient(const Data::Series &data, const Training::Task &task, size_t group_idx, Gradient &g) const;
    double class_likelihood_gradient(const Data::Set &data, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of the expected posterior probability of having at least one site with respect to transformed transition and emission probabilities */
    posterior_t posterior_gradient(const Data::Set &s, const Training::Task &task, size_t group_idx, matrix_t &transition_g, matrix_t &emission_g) const;
    posterior_gradient_t posterior_gradient(const Data::Seq &data, const Training::Task &task, size_t group_idx, matrix_t &transition_g, matrix_t &emission_g) const;


    /** (Log) likelihood gradient with respect to transformed transition probabilities */
    matrix_t transition_gradient(const matrix_t &T, const Training::Range &allowed) const;

    /** (Log) likelihood gradient with respect to transformed emission probabilities */
    matrix_t emission_gradient(const matrix_t &E, const Training::Range &allowed) const;


// -------------------------------------------------------------------------------------------
// Output methods
// -------------------------------------------------------------------------------------------
  public:
    std::string path2string(const StatePath &path) const;
    void to_dot(std::ostream &os, double minimum_transition=0.01) const;

// -------------------------------------------------------------------------------------------
// Consistency checking methods
// -------------------------------------------------------------------------------------------

    bool check_consistency_transitions() const;
    bool check_consistency_emissions() const;
    bool check_consistency() const;

    void switch_intermediate(bool new_state) { store_intermediate=new_state; };

    mask_t compute_mask(const Data::Collection &data) const;

// -------------------------------------------------------------------------------------------
// Methods for candidate generation in MCMC sampling
// -------------------------------------------------------------------------------------------
  protected:
    void swap_columns(std::mt19937 &rng);
    void modify_column(std::mt19937 &rng);
    void modify_transition(std::mt19937 &rng, double eps=1e-6);
    void add_columns(size_t n, std::mt19937 &rng);
    void add_column(size_t n, const std::vector<double> &e);
    void del_columns(size_t n, std::mt19937 &rng);
    void del_column(size_t n);

    struct History {
      public:
        size_t observation;
        size_t length;
        History() : observation(0), length(0) { };
    };

    void update_history(History &history, size_t obs) const;
    void update_history_front(History &history, size_t obs) const;
    size_t get_emission_index(size_t state, const History &history) const;

  public:
    bool is_motif_state(size_t state) const;
    bool is_motif_group(size_t group_idx) const;

    void register_dataset(const Data::Set &data, double class_prior, double motif_p1, double motif_p2);
  protected:
    Training::Range complementary_states(size_t group_idx) const;
};

std::ostream &operator<<(std::ostream& os, const HMM &hmm);

// NOTE: this assumes that the alphabet has size 4, which is okay for nucleic acids, but not for amino acids
inline size_t HMM::get_emission_index(size_t state, const History &history) const
{
  const size_t current_order_plus_1 = order[state] + 1;
  const size_t available_order_plus_1 = std::min<size_t>(history.length, current_order_plus_1);

  // truncate to available_order_plus_1 positions;
  // this is identity when available_order_plus_1 == order[state] + 1
  const size_t mask = (1 << (2 * available_order_plus_1)) - 1;
  const size_t obs = history.observation & mask;

  size_t res;
  if(available_order_plus_1 == current_order_plus_1)
    // return obs if the history size is limited by the order of the state 
    res = obs;
  else
    // otherwise add the corresponding order offset
    res = obs + order_offset[order[state] - available_order_plus_1];
  if(false) {
    std::cerr << "state = " << state << " history.observation = " << history.observation << " -> " << res << std::endl;
    std::cerr << "current_order_plus_1 = " << current_order_plus_1 << " available_order_plus_1 = " << available_order_plus_1 << " obs = " << obs << " mask = " << mask << " history.length = " << history.length << std::endl;
  }
  return(res);
}

confusion_matrix reduce(const vector_t &v, size_t group_idx, const Data::Series &data, bool word_stats);

#include "subhmm.hpp"

#endif

