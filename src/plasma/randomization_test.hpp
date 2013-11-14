
#include "data.hpp"
#include "stats.hpp"
#include "options.hpp"

bool randomization_test(const Seeding::DataCollection &collection, const Seeding::Stats::OccurrenceCounts &counts, size_t nr_tests, double score, const Seeding::Options &options, const Seeding::Objective &objective, size_t length, size_t degeneracy);

