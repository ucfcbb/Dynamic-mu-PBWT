#include <iostream>
#include <getopt.h>
#include <fstream>
#include <chrono>
#include <random>

#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
    std::cout << "Usage: dmupbwt [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --ref <path>\t Ref vcf file for panel" << std::endl;
    std::cout << "  -l, --load <path>\t Index file for panel" << std::endl;
    std::cout << "  -q, --insert <path>\t Vcf file of haplotypes to be inserted into the ref panel" << std::endl;
    std::cout << "  -o, --output log <path>\t txt file for logging the insertion times" << std::endl;
    std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

void Test_Insertion_EmptyPanel(string &ref_vcf_input, string &output_log, bool verbose) {
    DCPBWT dcpbwt;
    cout << "Testing Insertion...\n";
    vector<vector<bool>> alleles;
    clock_t START_query_read = clock();
    ReadQueryVCF(ref_vcf_input, alleles);
    //ReadQueryFile(ref_vcf_input.c_str(), alleles);
    auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
    cout << "Time to read haplotypes to be inserted: " << time_read_query << " s\n";

    // go through all query haplotypes
    ofstream out;
    out.open(output_log);
    out << "#No.\tHapID\tTime(ms)\tIndexMem(Bytes)\n";

    // randomize order of insertion
    std::random_device dev;
    std::mt19937 rng(dev());
    std::vector<int> insertion_order(alleles.size(), 0);
    iota(insertion_order.begin(), insertion_order.end(), 0);
    std::shuffle(insertion_order.begin(), insertion_order.end(), rng);

    clock_t START_INSERT_OVERALL = clock();
    for (unsigned int i = 0; i < insertion_order.size(); ++i) {
        auto begin = std::chrono::high_resolution_clock::now();
        dcpbwt.InsertSingleHaplotype(alleles[insertion_order[i]]);
        auto end = std::chrono::high_resolution_clock::now();
        auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        unsigned long long index_size = dcpbwt.get_memory_usage_bytes();
        out << i + 1 << "\t" << insertion_order[i] << "\t" << time_insert_per_hap << "\t" << index_size << "\n";
    }
    out.close();
    auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
    cout << "Inserted " << alleles.size() << " haplotypes.\n";
    cout << "Insertion took: " << time_insert << " s.\n";
    if (verbose)
        dcpbwt.PrintMemoryUsage(verbose);
    cout << "Avg. runs: " << dcpbwt.get_avg_runs() << "." << std::endl;
}

void Test_Insertion_RefPanel(string &ref_input, string &query_vcf_input, string &output_log, bool verbose, bool load) {
    if (load) {
        DCPBWT dcpbwt;
        std::ifstream in;
        clock_t START = clock();
        in.open(ref_input.c_str());
        dcpbwt.load(in);
        in.close();
        auto time_load = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Index loaded in " << time_load << " seconds." << std::endl;
        if (verbose)
            dcpbwt.PrintMemoryUsage(verbose);

        ofstream out;
        vector<double> insert_per_hap;
        out.open(output_log);

        cout << "Testing Insertion...\n";
        vector<vector<bool>> alleles;
        clock_t START_query_read = clock();
//  ReadQueryVCF(query_vcf_input, alleles);
        ReadQueryFile(query_vcf_input.c_str(), alleles);
        auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
        cout << "Time to read haplotypes to be inserted: " << time_read_query << " s\n";

        // go through all query haplotypes
        clock_t START_INSERT_OVERALL = clock();
        for (auto &allele : alleles) {
            auto begin = std::chrono::high_resolution_clock::now();
            dcpbwt.InsertSingleHaplotype(allele);
            auto end = std::chrono::high_resolution_clock::now();
            auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            out << time_insert_per_hap << "\n";
        }
        auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
        cout << "Inserted " << alleles.size() << " haplotypes.\n";
        cout << "Insertion took: " << time_insert << " s.\n";
        cout << "After insertion: \n";
        dcpbwt.PrintMemoryUsage(verbose);
        cout << "Avg. runs: " << dcpbwt.get_avg_runs() << "." << std::endl;
        out.close();
    } else {
        clock_t START = clock();
        DCPBWT dcpbwt(ref_input);
        auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
        cout << "Time to build: " << time_build << " secs.\n";
        if (verbose)
            dcpbwt.PrintMemoryUsage(verbose);

        ofstream out;
        vector<double> insert_per_hap;
        out.open(output_log);

        cout << "Testing Insertion...\n";
        vector<vector<bool>> alleles;
        clock_t START_query_read = clock();
//  ReadQueryVCF(query_vcf_input, alleles);
        ReadQueryFile(query_vcf_input.c_str(), alleles);
        auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
        cout << "Time to read haplotypes to be inserted: " << time_read_query << " s\n";

        // go through all query haplotypes
        clock_t START_INSERT_OVERALL = clock();
        for (auto &allele : alleles) {
            auto begin = std::chrono::high_resolution_clock::now();
            dcpbwt.InsertSingleHaplotype(allele);
            auto end = std::chrono::high_resolution_clock::now();
            auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            out << time_insert_per_hap << "\n";
        }
        auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
        cout << "Inserted " << alleles.size() << " haplotypes.\n";
        cout << "Insertion took: " << time_insert << " s.\n";
        dcpbwt.PrintMemoryUsage(verbose);
        cout << "Avg. runs: " << dcpbwt.get_avg_runs() << "." << std::endl;
        out.close();
    }
}

int main(int argc, char **argv) {
    if (argc == 1) {
        PrintHelp();
        exit(EXIT_SUCCESS);
    }
    std::string input_vcf, query_vcf, load_file;
    std::string output_log;
    bool load = false;
    bool verbose = false;

    int c = 0;
    while (true) {
        static struct option long_options[] = {
            {"ref", required_argument, nullptr, 'i'},
            {"load", required_argument, nullptr, 'l'},
            {"query", required_argument, nullptr, 'q'},
            {"output", required_argument, nullptr, 'o'},
            {"verbose", no_argument, nullptr, 'v'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:l:q:o:vh", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':input_vcf = optarg;
                break;
            case 'l':load_file = optarg;
                load = true;
                break;
            case 'q':query_vcf = optarg;
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

    if (query_vcf.empty()) {
        if (input_vcf.empty()) {
            cerr << "Input vcf isn't specified. Please specify this file with -i. " << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!filesystem::exists(input_vcf)) {
            cerr << "Input vcf : " << input_vcf << " doesn't exist!\n";
            exit(EXIT_FAILURE);
        }
        Test_Insertion_EmptyPanel(input_vcf, output_log, verbose);
    } else {
        if (load)
            Test_Insertion_RefPanel(load_file, query_vcf, output_log, verbose, load);
        else
            Test_Insertion_RefPanel(input_vcf, query_vcf, output_log, verbose, load);
    }
}
