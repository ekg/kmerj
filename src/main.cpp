#include <iostream>
#include <random>
#include <limits>
#include <cassert>
#include <chrono>
#include "args.hxx"
#include <endian.h>
#include "kmerj.hpp"

int main(int argc, char** argv) {

    args::ArgumentParser parser("kmerj (kmer histogram generator and comparator)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> in_file(parser, "FILE", "operate on this single input file", {'i', "in"});
    args::ValueFlag<std::string> in_file_a(parser, "FILE", "intersect kmers from this file", {'a', "in-a"});
    args::ValueFlag<std::string> in_file_b(parser, "FILE", "intersect kmers from this file", {'b', "in-b"});
    args::ValueFlag<std::string> db_file(parser, "FILE", "accumulate kmers in this temporary memory-mapped file", {'d', "db-file"});
    args::ValueFlag<uint64_t> kmer_size(parser, "N", "the number of bp in each kmer", {'k', "kmer-size"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

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

    uint64_t k = args::get(kmer_size);
    if (k > sizeof(uint64_t)*4) {
        std::cerr << "kmer size must be <= than " << sizeof(uint64_t)*4 << std::endl;
        exit(1);
    }

    if (!args::get(in_file).empty() && !args::get(db_file).empty()) {
        std::remove(args::get(db_file).c_str());
        mmmulti::set<uint64_t> ms(args::get(db_file));
        kmerj::kmerize(args::get(in_file), k, ms);
        uint64_t i = 0;
        uint64_t value_count = 0;
        uint64_t unique_value_count = 0;
        ms.for_each_value_count([&](const uint64_t& value, const uint64_t& count) {
                ++unique_value_count;
                value_count += count;
            });
        std::cerr << value_count << " values" << std::endl;
        std::cerr << unique_value_count << " uniques" << std::endl;
        std::remove(args::get(db_file).c_str());
    }

    if (!args::get(in_file_a).empty() && !args::get(in_file_b).empty()
        && !args::get(db_file).empty()) {
        std::string db_file_a = args::get(db_file) + ".a";
        std::string db_file_b = args::get(db_file) + ".b";
        std::remove(db_file_a.c_str());
        std::remove(db_file_b.c_str());
        mmmulti::set<uint64_t> multiset_a(db_file_a);
        mmmulti::set<uint64_t> multiset_b(db_file_b);
        kmerj::kmerize(args::get(in_file_a), k, multiset_a);
        kmerj::kmerize(args::get(in_file_b), k, multiset_b);
        auto handle_kmer = [&](const uint64_t& kmer,
                               const uint64_t& count_a,
                               const uint64_t& count_b) {
            std::cout << kmerj::unseq2bit(kmer, k) << "\t" << count_a << "\t" << count_b << std::endl;
        };
        kmerj::for_each_intersecting_kmer(multiset_a, multiset_b, handle_kmer);
        std::remove(db_file_a.c_str());
        std::remove(db_file_b.c_str());
    }

    return 0;

}
