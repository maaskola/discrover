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
 *       Filename:  data.hpp
 *
 *    Description:  Data structures to represent sets of FASTA sequences
 *
 *        Created:  Thu May 31 06:47:48 2012 +0200
 *
 *         Author:  Jonas Maaskola <jonas@maaskola.de>
 *
 * =====================================================================================
 */

#ifndef DATA_HPP
#define DATA_HPP

#include <map>
#include <unordered_map>
#include <unordered_set>
#include "specification.hpp"
#include "fasta.hpp"

std::string sha1hash(const std::string &s);

namespace Data {

  struct RemovalReport {
    RemovalReport(size_t n=0, size_t s=0);
    size_t nucleotides;
    size_t sequences;
  };

  RemovalReport operator+(const RemovalReport &a, const RemovalReport &b);
  RemovalReport &operator+=(RemovalReport &a, const RemovalReport &b);

  using mask_sub_t = std::unordered_map<std::string, std::vector<size_t>>;
  using mask_t = std::unordered_map<std::string, mask_sub_t>;

  namespace Basic {
    template <typename X>
      struct Set : public Specification::Set {

        // typedefs

        using seq_t = X;
        using seqs_t = std::vector<seq_t>;
        using iterator = typename seqs_t::iterator;
        using const_iterator = typename seqs_t::const_iterator;

        // constructors

        Set() :
          Specification::Set(),
          seq_size(0),
          set_size(0),
          sequences()
        { };
        Set(const Specification::Set &s, bool revcomp=false, size_t n_seq=0) :
          Specification::Set(s),
          seq_size(0),
          set_size(0),
          sequences()
        {
          read_fasta(path, sequences, revcomp, n_seq, is_shuffle);

          sha1 = compute_sha1();

          set_size = sequences.size();
          for(auto &seq: sequences)
            seq_size += seq.sequence.size();
        };
        template <typename Y>
          Set(const Set<Y> &set) :
            Specification::Set(set),
            seq_size(0),
            set_size(set.set_size),
            sequences(),
            sha1(set.sha1)
        {
          for(auto &seq: set) {
            seq_t s(seq);
            seq_size += s.sequence.size();
            sequences.push_back(s);
          }
        };

        std::string compute_sha1() const {
          std::string d;
          for(auto &s: sequences)
            d += s.sequence;
          return(sha1hash(d));
        }

        // member variables

        size_t seq_size, set_size;
        std::vector<seq_t> sequences;
        std::string sha1;

        RemovalReport mask(const mask_sub_t &mask) {
          RemovalReport report;
          for(auto &seq: sequences) {
            auto iter = mask.find(seq.definition);
            if(iter != end(mask)) {
              report.sequences++;
              report.nucleotides += seq.mask(iter->second);
            }
          }
          return(report);
        }

        RemovalReport drop_sequences(const mask_sub_t &mask) {
          const bool noisy_output = false;

          RemovalReport report;
          auto iter = sequences.rbegin();
          size_t idx = 0;
          while(iter != sequences.rend()) {
            if(noisy_output) {
              std::cerr << "idx = " << idx << std::endl;
              std::cerr << "Checking " << path << " " << iter->definition << " for dropping" << std::endl;
            }
            if(mask.find(iter->definition) != end(mask)) {
              if(noisy_output)
                std::cerr << "Dropping!" << std::endl;
              report.sequences++;
              report.nucleotides += iter->sequence.size();
              set_size--;
              seq_size -= iter->sequence.size(); // TODO: find out if this is correct for reverse complements
              auto i = iter;
              bool done = (++i) == sequences.rend();
              sequences.erase(--(iter++).base());
              if(done) break;
            } else
              iter++;
            if(noisy_output)
              std::cout << "idxB = " << idx++ << std::endl;
          }
          return(report);
        }
      };

    template <typename X>
      struct Contrast {

        // typedefs

        using set_t = X;
        using seq_t = typename set_t::seq_t;
        using sets_t = std::vector<set_t>;
        using iterator = typename sets_t::iterator;
        using const_iterator = typename sets_t::const_iterator;

        // constructors

        Contrast() : seq_size(0), set_size(0), sets(), name() { };
        Contrast(const std::string &name_, const Specification::Sets &specs, bool revcomp=false, size_t n_seq=0) :
          seq_size(0),
          set_size(0),
          sets(),
          name(name_)
        {
          for(auto &spec: specs)
            sets.push_back({spec, revcomp, n_seq});
          for(auto &s: sets) {
            seq_size += s.seq_size;
            set_size += s.set_size;
          }
        };
        template <typename Y>
          Contrast(const Contrast<Y> &contrast) :
            seq_size(0),
            set_size(contrast.set_size),
            sets(),
            name(contrast.name)
        {
          for(auto &x: contrast.sets) {
            set_t s(x);
            seq_size += s.seq_size;
            sets.push_back(s);
          }
        };

        // member variables

        size_t seq_size, set_size;
        sets_t sets;
        std::string name;

        // TODO add compression
        std::vector<std::string> save_shuffle_sequences(const std::string &stem, size_t &idx) const {
          std::vector<std::string> paths;
          for(auto dataset: sets)
            if(dataset.is_shuffle) {
              std::string path = stem + "_shuffle" + boost::lexical_cast<std::string>(idx++) + ".fa";
              paths.push_back(path);
              std::ofstream ofs(path.c_str());
              for(auto &seq: dataset)
                ofs << seq << std::endl;
            }
          return(paths);
        };

        RemovalReport mask(const Data::mask_t &mask) {
          RemovalReport report;
          for(auto &dataset: sets) {
            auto iter = mask.find(dataset.path);
            if(iter != end(mask))
              report += dataset.mask(iter->second);
          }
          return(report);
        }

        RemovalReport drop_sequences(const Data::mask_t &mask) {
          RemovalReport report;
          for(auto &dataset: sets) {
            auto iter = mask.find(dataset.path);
            if(iter != end(mask))
              report += dataset.drop_sequences(iter->second);
          }
          return(report);
        }
      };

    template <typename X>
      struct Collection {

        // typedefs

        using contrast_t = X;
        using set_t = typename contrast_t::set_t;
        using seq_t = typename set_t::seq_t;
        using contrasts_t = std::vector<contrast_t>;
        using iterator = typename contrasts_t::iterator;
        using const_iterator = typename contrasts_t::const_iterator;

        // constructors

        Collection() :
          seq_size(0),
          set_size(0),
          contrasts()
        { };
        Collection(const Specification::Sets &specs, bool revcomp=false, size_t n_seq=0) :
          seq_size(0),
          set_size(0),
          contrasts()
        {
          std::map<std::string,Specification::Sets> map;
          for(auto &spec: specs)
            map[spec.contrast].push_back(spec);
          for(auto &iter: map)
            contrasts.push_back(contrast_t(iter.first, iter.second, revcomp, n_seq));
          for(auto &s: contrasts) {
            seq_size += s.seq_size;
            set_size += s.set_size;
          }
        };
        template <typename Y>
          Collection(const Collection<Y> &coll) :
            seq_size(0),
            set_size(coll.set_size),
            contrasts()
        {
          for(auto &x: coll.contrasts) {
            contrast_t s(x);
            seq_size += s.seq_size;
            contrasts.push_back(s);
          }
        };

        // member variables

        size_t seq_size, set_size;
        contrasts_t contrasts;

        // methods

        iterator find(const std::string& s) {
          iterator iter = std::find_if(begin(contrasts), end(contrasts), [&s](const contrast_t &contrast) {
            return(contrast.name == s);
          });
          return(iter);
        }
        const_iterator find(const std::string& s) const {
          const_iterator iter = std::find_if(begin(contrasts), end(contrasts), [&s](const contrast_t &contrast) {
            return(contrast.name == s);
          });
          return(iter);
        }

        std::vector<std::string> save_shuffle_sequences(const std::string &path) const {
          std::vector<std::string> paths;
          size_t idx = 0;
          for(auto &contrast: contrasts)
            for(auto &x: contrast.save_shuffle_sequences(path, idx))
              paths.push_back(x);
          return(paths);
        };

        RemovalReport mask(const Data::mask_t &mask) {
          RemovalReport report;
          for(auto &contrast: contrasts)
            report += contrast.mask(mask);
          return(report);
        }

        RemovalReport drop_sequences(const Data::mask_t &mask) {
          RemovalReport report;
          for(auto &contrast: contrasts)
            report += contrast.drop_sequences(mask);
          return(report);
        }
      };

    template <typename X> typename Set<X>::iterator begin(Set<X> &set) { return(begin(set.sequences)); }
    template <typename X> typename Set<X>::iterator end(Set<X> &set) { return(end(set.sequences)); }
    template <typename X> typename Set<X>::const_iterator begin(const Set<X> &set) { return(begin(set.sequences)); }
    template <typename X> typename Set<X>::const_iterator end(const Set<X> &set) { return(end(set.sequences)); }

    template <typename X> typename Contrast<X>::iterator begin(Contrast<X> &contrast) { return(begin(contrast.sets)); }
    template <typename X> typename Contrast<X>::iterator end(Contrast<X> &contrast) { return(end(contrast.sets)); }
    template <typename X> typename Contrast<X>::const_iterator begin(const Contrast<X> &contrast) { return(begin(contrast.sets)); }
    template <typename X> typename Contrast<X>::const_iterator end(const Contrast<X> &contrast) { return(end(contrast.sets)); }

    template <typename X> typename Collection<X>::iterator begin(Collection<X> &collection) { return(begin(collection.contrasts)); }
    template <typename X> typename Collection<X>::iterator end(Collection<X> &collection) { return(end(collection.contrasts)); }
    template <typename X> typename Collection<X>::const_iterator begin(const Collection<X> &collection) { return(begin(collection.contrasts)); }
    template <typename X> typename Collection<X>::const_iterator end(const Collection<X> &collection) { return(end(collection.contrasts)); }
  }
}

namespace Seeding {

  struct Set : public Data::Basic::Set<Fasta::Entry> {
    Set();
    Set(const Specification::Set &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> Set(const Data::Basic::Set<X> &s) : Data::Basic::Set<Fasta::Entry>(s) { };
  };

  struct Contrast : public Data::Basic::Contrast<Set> {
    Contrast();
    Contrast(const std::string &name, const Specification::Sets &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> Contrast(const Data::Basic::Contrast<X> &s) : Data::Basic::Contrast<Set>(s) { };
  };

  struct Collection: public Data::Basic::Collection<Contrast> {
    Collection(const Specification::Sets &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> Collection(const Data::Basic::Collection<X> &c) : Data::Basic::Collection<Contrast>(c) { };
  };
};

#endif   /* ----- #ifndef DATA_HPP ----- */

