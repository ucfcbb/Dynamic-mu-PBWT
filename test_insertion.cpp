#include <iostream>
#include <getopt.h>
#include <fstream>
#include <chrono>

#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dcpbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
  std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
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
  auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
  cout << "Time to read query alleles: " << time_read_query << " s\n";

  // go through all query haplotypes
  ofstream out;
  out.open(output_log);
  vector<double> insert_per_hap;
  clock_t START_INSERT_OVERALL = clock();
  for (auto & allele : alleles) {
    auto begin = std::chrono::high_resolution_clock::now();
    dcpbwt.InsertSingleHaplotype(allele);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    out << time_insert_per_hap << "\n";
  }
  auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";
  dcpbwt.PrintMemoryUsage(verbose);
  out.close();
  cout << "Avg. runs: " << dcpbwt.get_avg_runs() << "." << std::endl;
}

void Test_Insertion_RefPanel(string &ref_vcf_input, string &query_vcf_input, string &output_log, bool verbose) {
  ofstream out;
  vector<double> insert_per_hap;
  out.open(output_log);
  clock_t START = clock();
  DCPBWT dcpbwt(ref_vcf_input, verbose);
  auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
  cout << "Time to build: " << time_build << " secs.\n";
  cout << "Testing Insertion...\n";
  vector<vector<bool>> alleles;
  clock_t START_query_read = clock();
  ReadQueryVCF(query_vcf_input, alleles);
  auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
  cout << "Time to read query alleles: " << time_read_query << " s\n";

  // go through all query haplotypes
  clock_t START_INSERT_OVERALL = clock();
  for (auto & allele : alleles) {
    auto begin = std::chrono::high_resolution_clock::now();
    dcpbwt.InsertSingleHaplotype(allele);
    auto end = std::chrono::high_resolution_clock::now(); 
    auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    out << time_insert_per_hap << "\n";
  }
  auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";
  dcpbwt.PrintMemoryUsage(verbose);
  cout << "Avg. runs: " << dcpbwt.get_avg_runs() << "." << std::endl;
  out.close();
}


int main(int argc, char **argv) {
  if (argc == 1) {
    PrintHelp();
    exit(EXIT_SUCCESS);
  }
  std::string input_vcf, query_vcf;
  std::string output_log;
  bool verbose = false;

  int c = 0;
  while (true) {
    static struct option long_options[] = {
      {"ref", required_argument, nullptr, 'i'},
      {"query", required_argument, nullptr, 'q'},
      {"output", required_argument, nullptr, 'o'},
      {"verbose", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "i:q:o:vh", long_options,
                    &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 'i':input_vcf = optarg;
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
  bool query = false;
  if (!filesystem::exists(input_vcf)) {
    cerr << "Input vcf : " << input_vcf << " doesn't exist!\n";
    exit(EXIT_FAILURE);
  }
  if (filesystem::exists(query_vcf)) {
    query = true;
  }

  if (query){
    Test_Insertion_RefPanel(input_vcf, query_vcf, output_log, verbose);
  } else {
    Test_Insertion_EmptyPanel(input_vcf, output_log, verbose);
  }
}
