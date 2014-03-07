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

// #define BOOST_UBLAS_MOVE_SEMANTICS // TODO find out if this helps

#include <boost/container/map.hpp>
#include <boost/container/flat_map.hpp>
#include <list>
#include <unordered_map>
#include <bitset>
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
  return(scalar_product(direction, gradient));
}

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
    HMM(Verbosity verbosity_, double pseudo_count_=1.0);
    HMM(const std::string &path, Verbosity verbosity_, double pseudo_count_=1.0);
    HMM(const HMM &hmm, bool copy_deep=true);

    typedef std::unordered_map<std::string, std::vector<size_t>> mask_sub_t;
    typedef std::unordered_map<std::string, mask_sub_t> mask_t;
    typedef boost::numeric::ublas::vector<size_t> StatePath;
    const static size_t max_motifs = 32;
    typedef std::bitset<max_motifs> bitmask_t;

  protected:
    /** The size of the alphabet. */
    const static size_t alphabet_size = 4;
    /** The number of emissions. */
    const static size_t n_emissions = alphabet_size;
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
      Specification::Set spec;
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

    void initialize_emissions();
    void initialize_transitions();
    void initialize_bg_transitions();
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
    bool check_enrichment(const Data::Contrast &contrast, const matrix_t &counts, size_t group_idx) const;

  public:
    /* Add a motif based on a IUPAC string. */
    size_t add_motif(const std::string &seq, double alpha, double exp_seq_len, double lambda, const std::string &name, const std::vector<size_t> &insertions, size_t pad_left, size_t pad_right);
    /* Add a motif based on a matrix. */
    size_t add_motif(const matrix_t &e, double exp_seq_len, double lambda, const std::string &name, std::vector<size_t> insertions, size_t pad_left, size_t pad_right);

    /* Add motifs of another HMM. */
    void add_motifs(const HMM &hmm, bool only_additional=false);

    /** Perform HMM training. */
    double train(const Data::Collection &col, const Training::Tasks &tasks, const hmm_options &options);
    /** Initialize HMM background with the Baum-Welch algorithm. */
    void initialize_bg_with_bw(const Data::Collection &col, const hmm_options &options);

  protected:
    double train_inner(const Data::Collection &col, const Training::Tasks &tasks, const hmm_options &options);

    /** Perform iterative HMM training, like re-estimation (expectation maximization) or gradient based. */
    void iterative_training(const Data::Collection &col, const Training::Tasks &tasks, const hmm_options &options);
    /** Perform one iteration of iterative HMM training. */
    bool perform_training_iteration(const Data::Collection &col, const Training::Tasks &tasks, const hmm_options &options, Training::State &ts);
    /** Perform one iteration of gradient training. */
    bool perform_training_iteration_gradient(const Data::Collection &col, const Training::Task &task, const hmm_options &options, int &center, double &score);
    /** Perform one iteration of re-estimation training. */
    bool perform_training_iteration_reestimation(const Data::Collection &col, const Training::Task &task, const hmm_options &options, double &score);
    /** Perform one iteration of class parameter re-estimation training. */
    bool reestimate_class_parameters(const Data::Collection &col, const Training::Task &task, const hmm_options &options, double &score);

    /** Perform HMM training for the background part of the model using Baum-Welch. */
    void train_background(const Data::Collection &col, const hmm_options &options);

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
  public:

    /** The emission information content of the HMM. */
    double information_content(size_t motif_idx) const;

    /** A routine that calculates the emission consensus. */
    std::string get_group_consensus(size_t idx, double threshold=0.1) const;

    /** A getter routine for the group names. */
    std::string get_group_name(size_t idx) const;

    /** A getter routine for the number of states. */
    size_t get_nstates() const;

    /** A getter routine for the number of groups, including the constitutive one. */
    size_t get_ngroups() const;

    /** A getter routine for the number of motifs, excluding the constitutive group. */
    size_t get_nmotifs() const;

    /** A getter routine for the length of a motif. */
    size_t get_motif_len(size_t motif_idx) const;

    /** A getter routine for the pseudo count. */
    double get_pseudo_count() const;

    size_t count_motif(const StatePath &path, size_t motif) const;

    /** Generate a random candidate variant for MCMC sampling */
    HMM random_variant(const hmm_options &options, std::mt19937 &rng) const;

    /** Prepare training tasks according to the objectives given in the options. */
    Training::Tasks define_training_tasks(const hmm_options &options) const;

  protected:

    size_t n_parameters() const;
    size_t non_zero_parameters(const Training::Targets &targets) const;

    friend std::ostream &operator<<(std::ostream& os, const HMM &hmm);

  public:
    double compute_score(const Data::Collection &col, const Measures::Continuous::Measure &measure, bool weighting, const std::vector<size_t> &present_motifs, const std::vector<size_t> &absent_motifs, const std::vector<size_t> &previous_motifs=std::vector<size_t>()) const;
    double compute_score_all_motifs(const Data::Collection &col, const Measures::Continuous::Measure &measure, bool weighting) const;

  protected:
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

  public:

    vector_t posterior_atleast_one(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;
    double viterbi(const Data::Seq &s, StatePath &path) const;
    posterior_t posterior_atleast_one(const Data::Seq &seq, bitmask_t present, bitmask_t absent) const;
    double expected_posterior(const Data::Seq &seq, bitmask_t present) const;

  protected:

    vector_t posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const;
    double   sum_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const;

    vector_t expected_posterior(const Data::Contrast &contrast, bitmask_t present) const;
    double   expected_posterior(const Data::Set &dataset, bitmask_t present) const;

  public:
    // posterior statistics of pairs of motifs
    struct pair_posterior_t {
      double log_likelihood;
      double posterior_first;
      double posterior_second;
      double posterior_both;
      double posterior_none;
    };

    typedef std::vector<pair_posterior_t> pair_posteriors_t;

  protected:

    pair_posteriors_t pair_posterior_atleast_one(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent, bitmask_t previous) const;

    pair_posteriors_t pair_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent, bitmask_t previous) const;
    pair_posterior_t sum_pair_posterior_atleast_one(const Data::Set &dataset, bitmask_t present, bitmask_t absent, bitmask_t previous) const;

// -------------------------------------------------------------------------------------------
// Generative and discriminative measures
// -------------------------------------------------------------------------------------------
// The logic for the present and absent arguments is as follows:
//  at least one of those specified as present has to be present
//  none of those specified as absent are to be present
// Furthermore, where present and absent groups are given by vectors, these are vectors of the indices of the respective groups.
// When they are to be given as masks, then these indices are transformed into bit masks.

  public:

    double log_likelihood(const Data::Collection &col) const;
    double log_likelihood(const Data::Contrast &contrast) const;
    double log_likelihood(const Data::Set &s) const;

    double class_likelihood(const Data::Contrast &contrast, const std::vector<size_t> &present_groups,  const std::vector<size_t> &absent_groups, bool compute_posterior) const;
    double class_likelihood(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent, bool compute_posterior) const;
    double class_likelihood(const Data::Set &dataset, bitmask_t present, bitmask_t absent, bool compute_posterior) const;

    // Discriminative measures, for groups specified by vectors of group indices
    /** The likelihood difference */
    double log_likelihood_difference(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups) const;
    /** The mutual information of condition and motif occurrence */
    double mutual_information(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups) const;
    /** The residual mutual information of condition and motif occurrence, as defined by Elemento et al in the method FIRE */
    double residual_mutual_information(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups, const std::vector<size_t> &previous_groups) const;
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups) const;
    /** The summed Matthew's correlation coefficients of all individual contrasts. */
    double matthews_correlation_coefficient(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups) const;
    /** The summed difference of expected occurrence frequencies. */
    double dips_tscore(const Data::Contrast &contrast, const std::vector<size_t> &present_groups) const;
    /** The summed difference of site occurrence. */
    double dips_sitescore(const Data::Contrast &contrast, const std::vector<size_t> &present_groups, const std::vector<size_t> &absent_groups) const;

  // protected:

    // Discriminative measures, for groups specified by bitmasks
    /** The likelihood difference */
    double log_likelihood_difference(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;
    /** The mutual information of condition and motif occurrence */
    double mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;
    /** The residual mutual information of condition and motif occurrence, as defined by Elemento et al in the method FIRE */
    double residual_mutual_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent, bitmask_t previous) const;
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;
    /** The summed Matthew's correlation coefficients of all individual contrasts. */
    double matthews_correlation_coefficient(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;
    /** The summed difference of expected occurrence frequencies. */
    double dips_tscore(const Data::Contrast &contrast, bitmask_t present) const;
    /** The summed difference of site occurrence. */
    double dips_sitescore(const Data::Contrast &contrast, bitmask_t present, bitmask_t absent) const;

    // Discriminative measures, for individual sets of sequences
    /** The mutual information of rank and motif occurrence. */
    double rank_information(const Data::Set &dataset, bitmask_t present, bitmask_t absent) const;

// -------------------------------------------------------------------------------------------
// Expectation-maximization type learning
// -------------------------------------------------------------------------------------------

  protected:

    /** Perform full expectation-maximization type learning */
    void reestimation(const Data::Collection &col, const Training::Task &task, const hmm_options &options);
    /** Perform one iteration of expectation-maximization type learning */
    double reestimationIteration(const Data::Collection &col, const Training::Task &task, const hmm_options &options);

    void M_Step(const matrix_t &T, const matrix_t &E, const Training::Targets &targets, const hmm_options &options);

    /** Perform one iteration of scaled Baum-Welch learning for a data collection */
    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Collection &col, const Training::Targets &targets, const hmm_options &options) const;
    /** Perform one iteration of scaled Baum-Welch learning for a data set */
    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Set &dataset, const Training::Targets &targets, const hmm_options &options) const;
    /** Perform one iteration of Viterbi learning for a data collection */
    double ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Collection &col, const Training::Targets &targets, const hmm_options &options);
    /** Perform one iteration of Viterbi learning for a data set */
    double ViterbiIteration(matrix_t &T, matrix_t &E, const Data::Set &dataset, const Training::Targets &targets, const hmm_options &options);

    double BaumWelchIteration(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const;
    double BaumWelchIteration_single(matrix_t &T, matrix_t &E, const Data::Seq &s, const Training::Targets &targets) const;

// -------------------------------------------------------------------------------------------
// Monte-Carlo Markov Chain inference
// -------------------------------------------------------------------------------------------
//
    /** Perform parallel tempering */
    std::vector<std::list<std::pair<HMM, double>>> mcmc(const Data::Collection &col, const Training::Task &task, const hmm_options &options);

// -------------------------------------------------------------------------------------------
// Gradient learning
// -------------------------------------------------------------------------------------------

    /** Perform gradient line searching with the desired objective function
     * using an approximate exponential step search.  **/
    std::pair<double, HMM> line_search(const Data::Collection &col, const Gradient &gradient, const Training::Task &task, int &center) const;
    /** Perform gradient line searching with the Moré-Thuente algorithm and the
     * desired objective function. */
    std::pair<double, HMM> line_search_more_thuente(const Data::Collection &col, const Gradient &gradient, double score, int &info, const Training::Task &task, const hmm_options &options) const;
    /** Auxiliary routine to build candidate HMM for a given direction and step
     * size. */
    HMM build_trial_model(const Gradient &gradient, double alpha, const Training::Task &task) const;

    /** Compute the gradient of the desired objective function */
    Gradient compute_gradient(const Data::Collection &col, double &score, const Training::Task &task, bool weighting) const;
    Gradient compute_gradient(const Data::Contrast &contrast, double &score, const Training::Task &task) const;
    double compute_gradient(const Data::Contrast &contrast, Gradient &gradient, const Training::Task &task, size_t group_idx) const;


// -------------------------------------------------------------------------------------------
// Forward and backward algorithms, and associated likelihood and posterior calculation
// -------------------------------------------------------------------------------------------

    /** The standard forward algorithm with scaling.
     *  The scaling vector is also determined. */
    matrix_t compute_forward_scaled(const Data::Seq &s, vector_t &scale) const;
    /** Computes only the scaling vector of the standard forward algorithm with scaling. */
    vector_t compute_forward_scale(const Data::Seq &s) const;

    /** The standard forward algorithm with pre-scaling.
     *  The scaling vector is assumed to be given. */
    matrix_t compute_forward_prescaled(const Data::Seq &s, const vector_t &scale) const;
    /** The standard backward algorithm with pre-scaling.
     *  The scaling vector is assumed to be given. */
    matrix_t compute_backward_prescaled(const Data::Seq &s, const vector_t &scale) const;

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
    double mutual_information_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of mutual information of rank and motif occurrence with respect to transformed transition and emission probabilities */
    double rank_information_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;
    double rank_information_gradient(const Data::Set &dataset, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of Matthew's correlation coefficient of class and motif occurrence with respect to transformed transition and emission probabilities */
    double matthews_correlation_coefficient_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of Chi-Square of class and motif occurrence with respect to transformed transition and emission probabilities */
    double chi_square_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of log likelihood difference with respect to transformed transition and emission probabilities */
    double log_likelihood_difference_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of site frequency difference with respect to transformed transition and emission probabilities */
    double site_frequency_difference_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of class likelihood with respect to transformed transition and emission probabilities */
    double class_likelihood_gradient(const Data::Contrast &contrast, const Training::Task &task, size_t group_idx, Gradient &g) const;
    double class_likelihood_gradient(const Data::Set &dataset, const Training::Task &task, size_t group_idx, Gradient &g) const;

    /** Gradient of the expected posterior probability of having at least one site with respect to transformed transition and emission probabilities */
    posterior_t posterior_gradient(const Data::Set &s, const Training::Task &task, size_t group_idx, matrix_t &transition_g, matrix_t &emission_g) const;
    posterior_gradient_t posterior_gradient(const Data::Seq &seq, const Training::Task &task, size_t group_idx, matrix_t &transition_g, matrix_t &emission_g) const;


    /** (Log) likelihood gradient with respect to transformed transition probabilities */
    matrix_t transition_gradient(const matrix_t &T, const Training::Range &allowed) const;

    /** (Log) likelihood gradient with respect to transformed emission probabilities */
    matrix_t emission_gradient(const matrix_t &E, const Training::Range &allowed) const;

  public:
    std::pair<HMM, std::map<size_t, size_t>> add_revcomp_motifs() const;

// -------------------------------------------------------------------------------------------
// Output methods
// -------------------------------------------------------------------------------------------
  public:
    std::string path2string_state(const StatePath &path) const;
    std::string path2string_group(const StatePath &path) const;
    void to_dot(std::ostream &os, double minimum_transition=0.01) const;

// -------------------------------------------------------------------------------------------
// Consistency checking methods
// -------------------------------------------------------------------------------------------

    bool check_consistency_transitions(double eps=1e-6) const;
    bool check_consistency_emissions(double eps=1e-6) const;
    bool check_consistency(double eps=1e-6) const;

    void switch_intermediate(bool new_state) { store_intermediate=new_state; };

    mask_t compute_mask(const Data::Collection &col) const;

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

    bool is_motif_state(size_t state) const;
  public:
    bool is_motif_group(size_t group_idx) const;

    void print_occurrence_table_header(std::ostream &out) const;
    void print_occurrence_table(const std::string &file, const Data::Seq &seq, const StatePath &path, std::ostream &out) const;

  protected:
    void register_dataset(const Data::Set &dataset, double class_prior, double motif_p1, double motif_p2);
    Training::Range complementary_states(size_t group_idx) const;
    Training::Range complementary_states_mask(bitmask_t present) const;
};


HMM::pair_posterior_t &operator+=(HMM::pair_posterior_t &one, const HMM::pair_posterior_t &two);
std::ostream &operator<<(std::ostream &os, const HMM::pair_posterior_t &one);
std::ostream &operator<<(std::ostream& os, const HMM &hmm);

#include "subhmm.hpp"

#endif

