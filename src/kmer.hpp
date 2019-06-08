#pragma once

#include <string>
#include <cassert>

namespace kmerj {

//bool has_non_ATGC(const std::string& s);
uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna);
std::string unseq2bit(uint64_t v, const uint64_t& size);

}
