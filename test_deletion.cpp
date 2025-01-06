#include <iostream>
#include <getopt.h>
#include <random>
#include <fstream>
#include <chrono>

#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
    std::cout << "Usage: dcpbwt [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input ref <path>\t input vcf file" << std::endl;
    std::cout << "  -l, --load <path>\t index file to be loaded" << std::endl;
    std::cout << "  -o, --output log <path>\t output file with deletion times" << std::endl;
    std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

void Test_Deletion(DCPBWT& dcpbwt, string &output_log, bool verbose) {
    ofstream out;
    out.open(output_log);
    out << "#No.\tHapID\tTime(ms)\tIndexMemBeforeDel(Bytes)\n";
    cout << "Testing Deletion...\n";
    std::random_device dev;
    std::mt19937 rng(dev());
    int cnt = 1;
    clock_t START_Del = clock();
    while (dcpbwt.M > 0) {
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0, dcpbwt.M - 1); // distribution in range [1, 6]
        unsigned int target_hap_id = dist6(rng);
        auto begin = std::chrono::high_resolution_clock::now();
        unsigned long long index_size = dcpbwt.get_memory_usage_bytes();
        dcpbwt.DeleteSingleHaplotype(target_hap_id);
        auto end = std::chrono::high_resolution_clock::now();
        auto del_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        out << cnt << "\t" << target_hap_id << "\t" << del_per_hap << "\t" << index_size << "\n";
        ++cnt;
    }
    dcpbwt.DeleteSingleHaplotype(0);
    auto time_del = (float) (clock() - START_Del) / CLOCKS_PER_SEC;
    cout << "Total Deletion took: " << time_del << " s.\n";
    out.close();
}

int main(int argc, char **argv) {
    if (argc == 1) {
        PrintHelp();
        exit(EXIT_SUCCESS);
    }
    std::string input_vcf, load_file;
    std::string output_log;
    bool verbose = false;

    int c = 0;
    while (true) {
        static struct option long_options[] = {
            {"ref", required_argument, nullptr, 'i'},
            {"load", required_argument, nullptr, 'l'},
            {"output", required_argument, nullptr, 'o'},
            {"verbose", no_argument, nullptr, 'v'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:l:o:vh", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':input_vcf = optarg;
                break;
            case 'l':load_file = optarg;
                break;
            case 'o':output_log = optarg;
                break;
            case 'v':verbose = true;
                break;
            case 'h':PrintHelp();
                exit(EXIT_SUCCESS);
            default:PrintHelp();
                exit(EXIT_FAILURE);
        }
    }
    if (input_vcf.empty() && load_file.empty()){
        cerr << "Input vcf or index file must be provided." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!load_file.empty()){
        if (!filesystem::exists(load_file)) {
            cerr << "Index file: " << load_file << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }
        DCPBWT dcpbwt;
        std::ifstream in;
        clock_t START = clock();
        in.open(load_file.c_str());
        dcpbwt.load(in);
        in.close();
        auto time_load = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Index loaded in " << time_load << " seconds." << std::endl;
        if (verbose){
            std::cout << "Total # of haplotypes = " << dcpbwt.M << std::endl;
            std::cout << "Total # of sites = " << dcpbwt.N << std::endl;
            dcpbwt.PrintMemoryUsage(verbose);
        }
        std::cout << "\n";
        Test_Deletion(dcpbwt, output_log, verbose);
        exit(EXIT_SUCCESS);
    }

    if (!input_vcf.empty()){
        if (!filesystem::exists(input_vcf)) {
            cerr << "Input vcf : " << input_vcf << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }
        clock_t START = clock();
        DCPBWT dcpbwt(input_vcf);
        auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
        cout << "Time to build: " << time_build << " secs.\n";
        if (verbose) {
            dcpbwt.PrintMemoryUsage(verbose);
        }
        std::cout << "\n";
        Test_Deletion(dcpbwt, output_log, verbose);
        exit(EXIT_SUCCESS);
    }
    exit(EXIT_SUCCESS);
}
