/* =====================================================================================
 * Copyright (c) 2012, Jonas Maaskola
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
 *       Filename:  plasma.hpp
 *
 *    Description:  Routines to find IUPAC regular expression type motifs
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef FIND_HPP
#define FIND_HPP

#include <future>
#include "options.hpp"
#include "plasma_stats.hpp"
#include "results.hpp"
#include "align.hpp"

namespace Seeding {

struct Plasma {
public:
  Options options;
  Collection collection;

private:
  bool needs_rebuilding;
  NucleotideIndex<size_t, size_t> index;

public:
  Plasma(const Options &options);
  Plasma(const Collection &collection_, const Options &opt);
  Results find_motifs(const Specification::Motif &motif,
                      const Objective &objective, bool doreport = true) const;
  void apply_mask(const Results &results);

private:
  Results find_seeds(size_t length, const Objective &objective,
                     Algorithm algorithm);
  Results find_plasma(size_t length, const Objective &objective,
                      size_t max_degeneracy,
                      const std::set<size_t> &degeneracies,
                      std::future<void> &index_rebuilt) const;
  Results find_external_dreme(size_t length, const Objective &objective,
                              size_t max_degeneracy,
                              const std::set<size_t> &degeneracies) const;
  Results find_mcmc(size_t length, const Objective &objective,
                    size_t max_degeneracy) const;
  rev_map_t determine_initial_candidates(
      size_t length, const Objective &objective, seq_type &best_motif,
      size_t &n_candidates, double &max_score, Results &results,
      const std::set<size_t> &degeneracies) const;
  Results find_all(const Specification::Motif &motif,
                   const Objective &objective) const;
  Results find_multiple(const Specification::Motif &motif,
                        const Objective &objective) const;
  void apply_mask(const std::string &motif);
  void apply_mask(const Result &result);
  std::future<void> rebuild_index();
};

void report(std::ostream &os, const Objective &objective,
            const std::string &motif, const Collection &collection,
            const Options &options);
void report(std::ostream &os, const Result &result,
            const Collection &collection, const Options &options);
void viterbi_dump(const std::string &motif, const Collection &collection,
                  std::ostream &out, const Options &options);

namespace Exception {
namespace Dreme {
struct OnlyBinaryContrast : std::exception {
  const char *what() const noexcept;
};
}
namespace Plasma {
struct NoObjectiveForMotif : std::exception {
  NoObjectiveForMotif(const std::string &token);
  const char *what() const noexcept;
  std::string token;
};
}
}
}

#endif /* ----- #ifndef FIND_HPP  ----- */
