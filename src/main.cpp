#include <iostream>
#include <random>
#include <limits>
#include <cassert>
#include <chrono>
#include "args.hxx"
#include <endian.h>
#include "gzstream.h"
#include "mmmultimap.hpp"
#include "mmmultiset.hpp"
#include "kmer.hpp"

int main(int argc, char** argv) {

    args::ArgumentParser parser("kmerj (kmer histogram generator and comparator)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> in_file(parser, "FILE", "take kmers from this input file", {'i', "in"});
    args::ValueFlag<std::string> db_file(parser, "FILE", "accumulate kmers in this temporary memory-mapped file", {'d', "db-file"});
    args::ValueFlag<uint64_t> kmer_size(parser, "N", "the number of bp in each kmer", {'k', "kmer-size"});
    //args::ValueFlag<uint64_t> max_val(parser, "N", "generate test data in the range [1,max_value]", {'M', "max-value"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    //args::ValueFlag<uint64_t> unique_value_tests(parser, "N", "number of unique value calls to make", {'u', "unique-vals"});
    //args::Flag test_multiset(parser, "multiset", "test the multiset", {'m', "test-multiset"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }

    if (!args::get(in_file).empty()) {
        std::remove(args::get(db_file).c_str());
        mmmulti::set<uint64_t> ms(args::get(db_file));
        igzstream in(args::get(in_file).c_str());
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
            assert(false);
        }
        uint64_t k = args::get(kmer_size);
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
            uint64_t end = seq.size()-k;
            const char* s = seq.c_str();
            // todo, make an omp task here
            for (uint64_t i = 0; i < end; ++i) {
                bool is_dna = true;
                uint64_t kint = kmerj::seq2bit(s+i, k, is_dna);
                if (is_dna) ms.append(kint);
            }
        }
        in.close();
        ms.index();
        uint64_t i = 0;
        uint64_t value_count = 0;
        uint64_t unique_value_count = 0;
        // exercise unique value search
        ms.for_each_value_count([&](const uint64_t& value, const uint64_t& count) {
                ++unique_value_count;
                value_count += count;
            });
        std::cerr << value_count << " values" << std::endl;
        std::cerr << unique_value_count << " uniques" << std::endl;
    }

    return 0;

}
