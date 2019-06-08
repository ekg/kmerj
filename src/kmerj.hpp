#pragma once

#include <string>
#include <cassert>
#include "gzstream.h"
#include "mmmultiset.hpp"

namespace kmerj {

uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna);
std::string unseq2bit(uint64_t v, const uint64_t& size);
double entropy(uint64_t v, const uint64_t& size);
double gc_rate(uint64_t v, const uint64_t& size);
void kmerize(const std::string& in_file, const uint64_t& k, mmmulti::set<uint64_t>& ms);
void for_each_intersecting_kmer(const mmmulti::set<uint64_t>& a, const mmmulti::set<uint64_t>& b,
                                const std::function<void(const uint64_t&, const uint64_t&, const uint64_t&)>& lambda);

}
