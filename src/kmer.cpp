#include "kmer.hpp"

namespace kmerj {

/*
bool has_non_ATGC(const std::string& s) {
    for (auto c : s) {
        switch (c) {
        case 'A': case 'a':
        case 'T': case 't':
        case 'G': case 'g':
        case 'C': case 'c':
            break;
        default:
            return false;
            break;
        }
    }
    return true;
}
*/

//uint64_t seq2bit(const std::string& s) {
uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna) {
    assert(k <= sizeof(uint64_t)*4);
    is_dna = true;
    uint64_t v = 0;
    for (uint64_t i = 0; i < k; ++i) {
        v <<= 2;
        switch (*(s+i)) {
        case 'A': case 'a':
            // 0b00, default case
            break;
        case 'T': case 't':
            v |= 0b01;
            break;
        case 'G': case 'g':
            v |= 0b10;
            break;
        case 'C': case 'c':
            v |= 0b11;
            break;
        default:
            is_dna = false;
            break;
        }
    }
    return v;
}

std::string unseq2bit(uint64_t v, const uint64_t& size) {
    assert(size <= sizeof(uint64_t)*4);
    std::string s(size, '\0');
    for (uint64_t i = size-1; i >= 0; --i) {
        // mask
        //uint64_t mask = 0b11 << i;
        switch (v & 0b11) {
        case 0b00:
            s[i] = 'A';
            break;
        case 0b01:
            s[i] = 'T';
            break;
        case 0b10:
            s[i] = 'G';
            break;
        case 0b11:
            s[i] = 'C';
            break;
        default:
            // impossible
            break;
        }
        v >>= 2;
    }
    return s;
}

}
