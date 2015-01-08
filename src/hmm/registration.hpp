/*
 * =====================================================================================
 *
 *       Filename:  registration.hpp
 *
 *    Description:  Store information about motif priors in classes
 *
 *        Created:  10/12/2014 08:33:02 AM
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef REGISTRATION_HPP
#define REGISTRATION_HPP

#include <unordered_map>
#include <string>
#include "basedefs.hpp"
#include "bitmask.hpp"
#include "../plasma/specification.hpp"
#include "../verbosity.hpp"

struct HMM;

struct Registration {
  struct Sample {
    Specification::Set spec;
    double class_prior;
    std::unordered_map<bitmask_t, double>
        motif_prior;  // the key is supposed to be the group index of the motif
    double get_motif_prior(bitmask_t present) const;
  };

  Registration(Verbosity verbosity = Verbosity::verbose);

  Verbosity verbosity;

  std::unordered_map<std::string, Sample>
      datasets;  // the keys are SHA1 hashes of the file contents

  double compute_marginal_motif_prior(bitmask_t present) const;
  double get_class_motif_prior(const std::string &sha1,
                               bitmask_t present) const;
  double get_class_prior(const std::string &sha1) const;

  void add_dataset(const Data::Set &dataset, double class_prior);
  void add_bitmask(const std::string name, bitmask_t present, double motif_p1,
                   double motif_p2);
};

#endif
