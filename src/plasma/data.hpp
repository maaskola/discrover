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
 *         Author:  Jonas Maaskola (JM), jonas@maaskola.de
 *
 * =====================================================================================
 */

#ifndef DATA_HPP
#define DATA_HPP

#include <map>
#include "specification.hpp"
#include "fasta.hpp"

namespace Data {
  namespace Basic {
    template <typename X>
      struct Set : public Specification::DataSet {

        // typedefs

        typedef X seq_t;
        typedef std::vector<seq_t> seqs_t;
        typedef typename seqs_t::iterator iterator;
        typedef typename seqs_t::const_iterator const_iterator;

        // constructors

        Set() :
          Specification::DataSet(),
          seq_size(0),
          set_size(0),
          sequences()
        { };
        Set(const Specification::DataSet &s, bool revcomp=false, size_t n_seq=0) :
          Specification::DataSet(s),
          seq_size(0),
          set_size(0),
          sequences()
        {
          read_fasta(path, sequences, revcomp, n_seq);
          set_size = sequences.size();
          for(auto &seq: sequences)
            seq_size += seq.sequence.size();
        };
        template <typename Y>
          Set(const Set<Y> &set) :
            Specification::DataSet(set),
            seq_size(0),
            set_size(set.set_size),
            sequences()
        {
          for(auto &seq: set) {
            seq_t s(seq);
            seq_size += s.sequence.size();
            sequences.push_back(s);
          }
        };

        // member variables

        size_t seq_size, set_size;
        std::vector<seq_t> sequences;
      };

    template <typename X>
      struct Series {

        // typedefs

        typedef X set_t;
        typedef typename set_t::seq_t seq_t;
        typedef std::vector<set_t> sets_t;
        typedef typename sets_t::iterator iterator;
        typedef typename sets_t::const_iterator const_iterator;

        // constructors

        Series() : seq_size(0), set_size(0), sets(), name() { };
        Series(const std::string &name_, const Specification::DataSets &specs, bool revcomp=false, size_t n_seq=0) :
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
          Series(const Series<Y> &series) :
            seq_size(0),
            set_size(series.set_size),
            sets(),
            name(series.name)
        {
          for(auto &x: series.sets) {
            set_t s(x);
            seq_size += s.seq_size;
            sets.push_back(s);
          }
        };

        // member variables

        size_t seq_size, set_size;
        sets_t sets;
        std::string name;
      };

    template <typename X>
      struct Collection {

        // typedefs

        typedef X series_t;
        typedef typename series_t::set_t set_t;
        typedef typename set_t::seq_t seq_t;
        typedef std::vector<series_t> serieses_t;
        typedef typename serieses_t::iterator iterator;
        typedef typename serieses_t::const_iterator const_iterator;

        // constructors

        Collection() :
          seq_size(0),
          set_size(0),
          series()
        { };
        Collection(const Specification::DataSets &specs, bool revcomp=false, size_t n_seq=0) :
          seq_size(0),
          set_size(0),
          series()
        {
          std::map<std::string,Specification::DataSets> map;
          for(auto &spec: specs)
            map[spec.series].push_back(spec);
          for(auto &iter: map)
            series.push_back(series_t(iter.first, iter.second, revcomp, n_seq));
          for(auto &s: series) {
            seq_size += s.seq_size;
            set_size += s.set_size;
          }
        };
        template <typename Y>
          Collection(const Collection<Y> &coll) :
            seq_size(0),
            set_size(coll.set_size),
            series()
        {
          for(auto &x: coll.series) {
            series_t s(x);
            seq_size += s.seq_size;
            series.push_back(s);
          }
        };

        // member variables

        size_t seq_size, set_size;
        serieses_t series;

        // methods

        iterator find(const std::string& s) {
          iterator iter = std::find_if(begin(series), end(series), [&s](const series_t &ser) {
            return(ser.name == s);
          });
          return(iter);
        }
        const_iterator find(const std::string& s) const {
          const_iterator iter = std::find_if(begin(series), end(series), [&s](const series_t &ser) {
            return(ser.name == s);
          });
          return(iter);
        }
      };

    template <typename X> typename Set<X>::iterator begin(Set<X> &set) { return(begin(set.sequences)); }
    template <typename X> typename Set<X>::iterator end(Set<X> &set) { return(end(set.sequences)); }
    template <typename X> typename Set<X>::const_iterator begin(const Set<X> &set) { return(begin(set.sequences)); }
    template <typename X> typename Set<X>::const_iterator end(const Set<X> &set) { return(end(set.sequences)); }

    template <typename X> typename Series<X>::iterator begin(Series<X> &series) { return(begin(series.sets)); }
    template <typename X> typename Series<X>::iterator end(Series<X> &series) { return(end(series.sets)); }
    template <typename X> typename Series<X>::const_iterator begin(const Series<X> &series) { return(begin(series.sets)); }
    template <typename X> typename Series<X>::const_iterator end(const Series<X> &series) { return(end(series.sets)); }

    template <typename X> typename Collection<X>::iterator begin(Collection<X> &collection) { return(begin(collection.series)); }
    template <typename X> typename Collection<X>::iterator end(Collection<X> &collection) { return(end(collection.series)); }
    template <typename X> typename Collection<X>::const_iterator begin(const Collection<X> &collection) { return(begin(collection.series)); }
    template <typename X> typename Collection<X>::const_iterator end(const Collection<X> &collection) { return(end(collection.series)); }
  }
}

#include <unordered_set>
#include <unordered_map>

namespace Plasma {
  typedef std::unordered_map<std::string, std::vector<size_t>> mask_sub_t;
  typedef std::unordered_map<std::string, mask_sub_t> mask_t;

  struct RemovalReport {
    RemovalReport(size_t n=0, size_t s=0);
    size_t nucleotides;
    size_t sequences;
  };

  RemovalReport operator+(const RemovalReport &a, const RemovalReport &b);
  RemovalReport &operator+=(RemovalReport &a, const RemovalReport &b);

  struct DataSet : public Data::Basic::Set<Fasta::Entry> {
    DataSet();
    DataSet(const Specification::DataSet &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> DataSet(const Data::Basic::Set<X> &s) : Data::Basic::Set<Fasta::Entry>(s) { };
    RemovalReport mask(const mask_sub_t &mask, char mask_symbol='n');
    RemovalReport drop_sequences(const std::unordered_set<std::string> &ids);
  };

  struct DataSeries : public Data::Basic::Series<DataSet> {
    DataSeries();
    DataSeries(const std::string &name, const Specification::DataSets &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> DataSeries(const Data::Basic::Series<X> &s) : Data::Basic::Series<DataSet>(s) { };
    RemovalReport mask(const mask_t &mask);
    RemovalReport drop_sequences(std::map<std::string,std::unordered_set<std::string>> &ids);
  };

  struct DataCollection: public Data::Basic::Collection<DataSeries> {
    DataCollection(const Specification::DataSets &s, bool revcomp=false, size_t n_seq=0);
    template <typename X> DataCollection(const Data::Basic::Collection<X> &c) : Data::Basic::Collection<DataSeries>(c) { };
    RemovalReport mask(const mask_t &mask);
    RemovalReport drop_sequences(std::map<std::string,std::unordered_set<std::string>> &ids);
  };
};

#endif   /* ----- #ifndef DATA_HPP ----- */

