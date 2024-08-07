#include <iostream>
#include <getopt.h>

#include "dcpbwt.h"


void PrintHelp() {
    std::cout << "Usage: dcpbwt [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input_file <path>\t vcf file for panel" << std::endl;
    std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

int main(int argc, char** argv){
    if (argc == 1) {
        PrintHelp();
        exit(EXIT_SUCCESS);
    }
    std::string ref_vcf_input;
    std::string output;
    bool verbose = false;

    int c = 0;
    while (true){
        static struct option long_options[] = {
                {"input",   required_argument, nullptr, 'i'},
                {"verbose",   no_argument, nullptr, 'v'},
                {"help",    no_argument,       nullptr, 'h'},
                {nullptr, 0,                   nullptr, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:vh", long_options,
                        &option_index);

        if (c == -1) {
                break;
        }

        switch (c) {
            case 'i':
                ref_vcf_input = optarg;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                PrintHelp();
                exit(EXIT_SUCCESS);
            default:
                PrintHelp();
                exit(EXIT_FAILURE);
        }
    }

    DCPBWT dcpbwt(ref_vcf_input, verbose);
    cout << dcpbwt.combined[0].size() << "\n";
    cout << dcpbwt.run_head(0, 0) << "\n";
    cout << dcpbwt.run_head(0, 1) << "\n";
    cout << dcpbwt.run_head(0, 2) << "\n";
    cout << dcpbwt.run_head(0, 3) << "\n";
    return (EXIT_SUCCESS);
}
