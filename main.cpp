#include <iostream>
#include <getopt.h>
#include <random>
#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dmupbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
  std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
  std::cout << "  -o, --output file <path>\t output file of matches" << std::endl;
  std::cout << "  -l, --input length <int>\t of match (in sites)" << std::endl;
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
  std::string output_file;
  unsigned int length = 0;
  bool verbose = false;

  int c = 0;
  while (true) {
    static struct option long_options[] = {
      {"reference", required_argument, nullptr, 'i'},
      {"query", required_argument, nullptr, 'q'},
      {"length", required_argument, nullptr, 'l'},
      {"output", required_argument, nullptr, 'o'},
      {"verbose", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "i:q:l:o:bvh", long_options,
                    &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 'i':ref_vcf_input = optarg;
        break;
      case 'q':query_vcf_input = optarg;
        break;
      case 'l':length = std::stoi(optarg);
        break;
      case 'o':output_file = optarg;
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
  if (ref_vcf_input.empty()) {
    cerr << "Input ref vcf file is not provided! Must provide this file.\n";
    exit(EXIT_FAILURE);
  }

  if (!query_vcf_input.empty()) {
    query = true;
  } else {
      clock_t START = clock();
      DCPBWT dcpbwt(ref_vcf_input);
      auto time_build = (float)(clock() - START)/CLOCKS_PER_SEC;
      std::cout << "Build took " << time_build << " seconds." << std::endl;
      if (verbose){
        dcpbwt.PrintMemoryUsage(verbose);
      }
      return EXIT_SUCCESS;
  }

  // Long Match Query
  if (query && length > 0) {
    if (output_file.empty()) {
      cerr << "Output file to store matches missing! Please specify an output file.\n";
      exit(EXIT_FAILURE);
    }
    if (!filesystem::exists(ref_vcf_input)) {
      cerr << "Specified ref vcf : " << ref_vcf_input << " doesn't exist!\n";
      exit(EXIT_FAILURE);
    }
    if (!filesystem::exists(query_vcf_input)) {
      cerr << "Specified query vcf : " << query_vcf_input << " doesn't exist!\n";
      exit(EXIT_FAILURE);
    }
    DCPBWT dcpbwt(ref_vcf_input);
    if (verbose){
      dcpbwt.PrintMemoryUsage(verbose);
    }
    dcpbwt.long_match_query(query_vcf_input, output_file, static_cast<unsigned int>(length), verbose);
  }
  return EXIT_SUCCESS;
}

