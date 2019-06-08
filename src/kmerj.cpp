#include "kmerj.hpp"

namespace kmerj {

uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna) {
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
    for (int64_t i = size-1; i >= 0; --i) {
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
            break;
        }
        v >>= 2;
    }
    return s;
}

void kmerize(const std::string& in_file, const uint64_t& k, mmmulti::set<uint64_t>& ms) {
    igzstream in(in_file.c_str());
    bool input_is_fasta=false, input_is_fastq=false;
    // look at the first character to determine if it's fastq or fasta
    std::string line; // line buffer
    std::getline(in, line);
    if (line[0] == '>') {
        input_is_fasta = true;
    } else if (line[0] == '@') {
        input_is_fastq = true;
    } else {
        std::cerr << "unknown file format" << std::endl;
        exit(1);
    }
#pragma omp parallel default(none) shared(in, ms, input_is_fasta, input_is_fastq, k, line)
#pragma omp single
    while (in.good()) {
        std::string seq;
        // get the sequence
        if (input_is_fasta) {
            while (std::getline(in, line)) {
                if (line[0] == '>') {
                    // this is the header of the next sequence
                    break;
                } else {
                    seq.append(line);
                }
            }
        } else if (input_is_fastq) {
            std::getline(in, seq); // sequence
            std::getline(in, line); // delimiter
            std::getline(in, line); // quality
            std::getline(in, line);
        }
        // add each kmer in the seq to our set
#pragma omp task default(none) firstprivate(seq) shared(ms, k)
        {
            uint64_t end = seq.size()-k;
            const char* s = seq.c_str();
            // todo, make an omp task here
            for (uint64_t i = 0; i < end; ++i) {
                bool is_dna = true;
                uint64_t kint = seq2bit(s+i, k, is_dna);
                if (is_dna) ms.append(kint);
                /*
                  if (is_dna && unseq2bit(kint, k) != seq.substr(i, k)) {
                  std::cerr << "for " << kint << " " << unseq2bit(kint, k) << " != " << seq.substr(i, k) << std::endl;
                  exit(1);
                  }
                */
            }
        }
    }
    in.close();
    ms.index();
}

void for_each_intersecting_kmer(const mmmulti::set<uint64_t>& a, const mmmulti::set<uint64_t>& b,
                                const std::function<void(const uint64_t&, const uint64_t&, const uint64_t&)>& lambda) {
    uint64_t a_idx = 0;
    uint64_t b_idx = 0;
    uint64_t a_size = a.size();
    uint64_t b_size = b.size();
    if (!a_size || !b_size) return;
    uint64_t a_kmer = a.read_value(a_idx++);
    uint64_t b_kmer = b.read_value(b_idx++);
    auto count_common_kmers = [&](void) {
        while (a_idx < a_size && b_idx < b_size
               && a_kmer == b_kmer) {
            uint64_t common_kmer = a_kmer;
            uint64_t a_start = a_idx;
            uint64_t b_start = b_idx;
            while (a_idx < a_size && a_kmer == common_kmer) {
                a_kmer = a.read_value(a_idx++);
            }
            while (b_idx < b_size && b_kmer == common_kmer) {
                b_kmer = b.read_value(b_idx++);
            }
            uint64_t a_count = a_idx - a_start;
            uint64_t b_count = b_idx - b_start;
            lambda(common_kmer, a_count, b_count);
        }
    };
    while (a_idx < a_size && b_idx < b_size) {
        while (a_idx < a_size && a_kmer < b_kmer) {
            a_kmer = a.read_value(a_idx++);
        }
        count_common_kmers();
        while (b_idx < b_size && b_kmer < a_kmer) {
            b_kmer = b.read_value(b_idx++);
        }
        count_common_kmers();
    }
}

}
