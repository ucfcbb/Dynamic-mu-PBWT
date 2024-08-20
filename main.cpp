#include <iostream>
#include <getopt.h>

#include "dcpbwt.h"
#include <random>

#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dcpbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input_file <path>\t vcf file for panel" << std::endl;
  std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
  std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

void Test_UV(DCPBWT &dcpbwt) {
  unsigned int query_idx = 0;
  auto uv = dcpbwt.uv_trick(0, query_idx);
  cout << "For i = " << query_idx << ": u = " << uv.first << ", v = " << uv.second << "\n";

  query_idx = 5;
  uv = dcpbwt.uv_trick(0, query_idx);
  cout << "For i = " << query_idx << ": u = " << uv.first << ", v = " << uv.second << "\n";

  query_idx = 6;
  uv = dcpbwt.uv_trick(0, query_idx);
  cout << "For i = " << query_idx << ": u = " << uv.first << ", v = " << uv.second << "\n";

  query_idx = 10;
  uv = dcpbwt.uv_trick(0, query_idx);
  cout << "For i = " << query_idx << ": u = " << uv.first << ", v = " << uv.second << "\n";

  query_idx = 17;
  uv = dcpbwt.uv_trick(0, query_idx);
  cout << "For i = " << query_idx << ": u = " << uv.first << ", v = " << uv.second << "\n";
}

void Test_BottomUp_Delete(DCPBWT &dcpbwt){
  cout << "Testing Deletion BottomUp i.e. from last haplotype to first haplotype in order...\n";
  for(int i = dcpbwt.M - 1; i >= 0; --i){
    dcpbwt.DeleteSingleHaplotype(i);
    cout << "Deleted hap " << i << "\n";
  }
  dcpbwt.DeleteSingleHaplotype(0);
}

void Test_TopDown_Delete(DCPBWT &dcpbwt){
  cout << "Testing Deletion TopDown i.e. from first haplotype to last haplotype in order...\n";
  for(unsigned int i = 0; i < dcpbwt.M; ++i){
    dcpbwt.DeleteSingleHaplotype(i);
    cout << "Deleted hap " << i << "\n";
  }
  dcpbwt.DeleteSingleHaplotype(0);
}

void Test_RandomDelete(DCPBWT& dcpbwt){
  cout << "Testing Random Deletion ...\n";
  std::mt19937 gen( std::random_device{}() );
  while (!dcpbwt.haplotype_ids.empty()){
    unsigned int element = 0;
    std::sample( dcpbwt.haplotype_ids.begin(), dcpbwt.haplotype_ids.end(), &element, 1, gen );
    cout << element << "\n";
    if (dcpbwt.haplotype_ids.find(element) != dcpbwt.haplotype_ids.end()){
      dcpbwt.DeleteSingleHaplotype(element);
      cout << "Deleted hap " << element << "\n";
    }
  }
//  dcpbwt.DeleteSingleHaplotype(0);

}

void Test_Insertion(string& ref_vcf_input, string& query_vcf_input, bool verbose){
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
  clock_t START_INSERT = clock();
  for(int i = 0; i < alleles.size(); ++i) {
    dcpbwt.InsertSinglelHaplotype(alleles[i]);
    cout << "Inserted query hap: " << i << "\n";
  }
  auto time_insert = (float) (clock() - START_INSERT) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";
}

int main(int argc, char **argv) {
  if (argc == 1) {
    PrintHelp();
    exit(EXIT_SUCCESS);
  }
  std::string ref_vcf_input;
  std::string query_vcf_input;
  std::string output;
  bool verbose = false;

  int c = 0;
  while (true) {
    static struct option long_options[] = {
      {"reference", required_argument, nullptr, 'i'},
      {"query", required_argument, nullptr, 'q'},
      {"verbose", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "i:q:vh", long_options,
                    &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 'i':ref_vcf_input = optarg;
        break;
      case 'q':query_vcf_input = optarg;
        break;
      case 'v':verbose = true;
        break;
      case 'h':PrintHelp();
        exit(EXIT_SUCCESS);
      default:PrintHelp();
        exit(EXIT_FAILURE);
    }
  }

  // Test_UV(dcpbwt);

//  vector<bool> query{false, true, false, false, true, false, true, false, false, false, true, true, true, false, true};
//  dcpbwt.InsertSinglelHaplotype(query);

//  Test_BottomUp_Delete(dcpbwt);
//  Test_TopDown_Delete(dcpbwt);
//  Test_RandomDelete(dcpbwt);
  Test_Insertion(ref_vcf_input, query_vcf_input, verbose);


  // Test_Insert();
  return (EXIT_SUCCESS);
}
