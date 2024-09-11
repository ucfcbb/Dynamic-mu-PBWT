#include <iostream>
#include <getopt.h>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dcpbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
  std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
  std::cout << "  -o, --output log <path>\t txt file for log" << std::endl;
  std::cout << "  -l, --input length <int>\t of match" << std::endl;
  std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
  std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

void Test_Insertion(string &ref_vcf_input, string &query_vcf_input, string &output_log, bool verbose) {
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
  for (int i = 0; i < alleles.size(); ++i) {
    auto begin = std::chrono::high_resolution_clock::now();
    dcpbwt.InsertSingleHaplotype(alleles[i]);
    auto end = std::chrono::high_resolution_clock::now(); 
    auto time_insert_per_hap = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    //insert_per_hap.push_back(time_insert_per_hap);
    out << time_insert_per_hap << "\n";
  }
  auto time_insert = (float) (clock() - START_INSERT_OVERALL) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";

  //if (out.is_open()){
  //  for(auto val: insert_per_hap){
  //    out << std::setprecision(10) << val << "\n";
  //  }
  //}
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
  if (!filesystem::exists(input_vcf)) {
    cerr << "Input vcf : " << input_vcf << " doesn't exist!\n";
    exit(EXIT_FAILURE);
  }
  if (!filesystem::exists(query_vcf)) {
    cerr << "Query vcf : " << query_vcf << " doesn't exist!\n";
    exit(EXIT_FAILURE);
  }

  Test_Insertion(input_vcf, query_vcf, output_log,verbose);
}
