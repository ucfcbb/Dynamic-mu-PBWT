#include <iostream>
#include <getopt.h>
#include <fstream>
#include <random>
#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
    std::cout << "Usage: dmupbwt [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
    std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
    std::cout << "  -o, --output file <path>\t output file of matches" << std::endl;
    std::cout << "  -L, --input length <int>\t of match (in sites)" << std::endl;
    std::cout << "  -l, --load index " << std::endl;
    std::cout << "  -s, --save <path>\t save the built dynamic mu-pbwt index" << std::endl;
    std::cout << "  -v, --verbose <path>\t show detail information and memory usage" << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        PrintHelp();
        exit(EXIT_SUCCESS);
    }
    std::string ref_vcf_input;
    std::string query_vcf_input;
    std::string load_file;
    std::string output_file;
    std::string save_index_file;
    unsigned int length = 0;
    bool verbose = false;

    int c = 0;
    while (true) {
        static struct option long_options[] = {
            {"reference", required_argument, nullptr, 'i'},
            {"query", required_argument, nullptr, 'q'},
            {"length", required_argument, nullptr, 'L'},
            {"output", required_argument, nullptr, 'o'},
            {"load", required_argument, nullptr, 'l'},
            {"save", required_argument, nullptr, 's'},
            {"verbose", no_argument, nullptr, 'v'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:q:L:o:s:l:vh", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':ref_vcf_input = optarg;
                break;
            case 'q':query_vcf_input = optarg;
                break;
            case 'L':length = std::stoi(optarg);
                break;
            case 'o':output_file = optarg;
                break;
            case 's':save_index_file = optarg;
                break;
            case 'l':load_file = optarg;
                break;
            case 'v':verbose = true;
                break;
            case 'h':PrintHelp();
                exit(EXIT_SUCCESS);
            default:PrintHelp();
                exit(EXIT_FAILURE);
        }
    }

    bool query = false;
    if (ref_vcf_input.empty() && load_file.empty()) {
        cerr << "Input ref vcf file (or index file) is not provided! Must provide this file.\n";
        exit(EXIT_FAILURE);
    }

    if (!query_vcf_input.empty()) {
        if (!filesystem::exists(query_vcf_input)) {
            cerr << "Specified query vcf : " << query_vcf_input << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }
        if (length <= 0){
            cerr << "Please specify length > 0.\n";
            exit(EXIT_FAILURE);
        }
        if (output_file.empty()) {
            cerr << "Output file to store matches missing! Please specify an output file to store matches.\n";
            exit(EXIT_FAILURE);
        }
        query = true;
    }

    if (!load_file.empty()) {
        if (!filesystem::exists(load_file)) {
            cerr << "Specified index file : " << load_file << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }

        DCPBWT dcpbwt;
        std::ifstream load;
        clock_t START = clock();
        load.open(load_file.c_str());
        dcpbwt.load(load);
        load.close();
        auto time_load = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Index loaded in " << time_load << " seconds." << std::endl;
        if (verbose) {
            dcpbwt.PrintMemoryUsage(verbose);
        }

        // Long Match Query
        if (query){
            dcpbwt.long_match_query(query_vcf_input, output_file, static_cast<unsigned int>(length), verbose);
        }
    } else {
        if (!filesystem::exists(ref_vcf_input)) {
            cerr << "Specified ref vcf : " << ref_vcf_input << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }

        clock_t START = clock();
        DCPBWT dcpbwt(ref_vcf_input);
        auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Build took " << time_build << " seconds." << std::endl;

        if (!save_index_file.empty()){
            std::ofstream outstream;
            outstream.open(save_index_file.c_str());
            if (!outstream.is_open()) {
                std::cerr << "Error opening file!" << std::endl;
                return 1;
            }
            auto written_bytes = dcpbwt.serialize(outstream);
            outstream.close();
            std::cout << "Wrote " << written_bytes << " bytes of index to " << save_index_file << "." << std::endl;
        }
        if (verbose) {
            dcpbwt.PrintMemoryUsage(verbose);
        }

        // Long Match Query
        if (query){
            dcpbwt.long_match_query(query_vcf_input, output_file, static_cast<unsigned int>(length), verbose);
        }
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}

