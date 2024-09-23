#include <iostream>
#include <getopt.h>
#include <random>

#include "dcpbwt.h"
#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dcpbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
  std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
  std::cout << "  -o, --output file <path>\t output file of matches" << std::endl;
  std::cout << "  -l, --input length <int>\t of match" << std::endl;
  std::cout << "  -v, --verbose <path>\t show detail information" << std::endl;
  std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

void Check_Insertion(DCPBWT &dcpbwt_inserted, DCPBWT &dcpbwt_build) {
  assert(dcpbwt_inserted.N == dcpbwt_build.N);
  assert(dcpbwt_inserted.M == dcpbwt_build.M);

  // Test if size of each of the spsis match
  for (unsigned int col = 0; col <= dcpbwt_build.N; ++col) {
    assert(dcpbwt_inserted.columns[col].num_zeros == dcpbwt_build.columns[col].num_zeros);
    assert(dcpbwt_inserted.columns[col].start_with_zero == dcpbwt_build.columns[col].start_with_zero);
    assert(dcpbwt_inserted.columns[col].ones.size() == dcpbwt_build.columns[col].ones.size());
    assert(dcpbwt_inserted.columns[col].zeros.size() == dcpbwt_build.columns[col].zeros.size());
    assert(dcpbwt_inserted.columns[col].combined.size() == dcpbwt_build.columns[col].combined.size());
    assert(dcpbwt_inserted.columns[col].pref_samples_beg.size() == dcpbwt_build.columns[col].pref_samples_beg.size());
    assert(dcpbwt_inserted.columns[col].pref_samples_end.size() == dcpbwt_build.columns[col].pref_samples_end.size());
    assert(dcpbwt_inserted.columns[col].div_samples_beg.size() == dcpbwt_build.columns[col].div_samples_beg.size());
  }

  for (unsigned int col = 0; col <= dcpbwt_build.N; ++col) {
    // Test ones
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].ones.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].ones.at(i) == dcpbwt_build.columns[col].ones.at(i));
    }
    std::cout << "Passed ones in col " << col << std::endl;
    // Test zeros
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].zeros.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].zeros.at(i) == dcpbwt_build.columns[col].zeros.at(i));
    }
    std::cout << "Passed zeros in col " << col << std::endl;
    // Test combined
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].combined.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].combined.at(i) == dcpbwt_build.columns[col].combined.at(i));
    }
    std::cout << "Passed combined in col " << col << std::endl;
    // Test pref_sample_beg
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].pref_samples_beg.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].pref_samples_beg.at(i) == dcpbwt_build.columns[col].pref_samples_beg.at(i));
    }
    std::cout << "Passed pref_samples_beg in col " << col << std::endl;
    // Test pref_sample_end
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].pref_samples_end.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].pref_samples_end.at(i) == dcpbwt_build.columns[col].pref_samples_end.at(i));
    }
    std::cout << "Passed pref_samples_end in col " << col << std::endl;
    // Test div_sample_beg
    for (unsigned int i = 0; i < dcpbwt_build.columns[col].div_samples_beg.size(); ++i) {
      assert(dcpbwt_inserted.columns[col].div_samples_beg.at(i) == dcpbwt_build.columns[col].div_samples_beg.at(i));
    }
    std::cout << "Passed div_samples_beg in col " << col << std::endl;
  }

  std::cout << "Test on phi structures" << std::endl;
  for (unsigned int hap = 0; hap < dcpbwt_build.M; ++hap) {
    assert(dcpbwt_inserted.phi->phi_vec[hap].size() == dcpbwt_build.phi->phi_vec[hap].size());
    assert(dcpbwt_inserted.phi->phi_inv_vec[hap].size() == dcpbwt_build.phi->phi_inv_vec[hap].size());
    assert(dcpbwt_inserted.phi->phi_supp[hap].size() == dcpbwt_build.phi->phi_supp[hap].size());
    assert(dcpbwt_inserted.phi->phi_inv_supp[hap].size() == dcpbwt_build.phi->phi_inv_supp[hap].size());
    assert(dcpbwt_inserted.phi->phi_supp_lcp[hap].size() == dcpbwt_build.phi->phi_supp_lcp[hap].size());

    for (unsigned int j = 0; j < dcpbwt_build.phi->phi_vec[hap].size(); ++j) {
      assert(dcpbwt_inserted.phi->phi_vec[hap].at(j) == dcpbwt_build.phi->phi_vec[hap].at(j));
      assert(dcpbwt_inserted.phi->phi_inv_vec[hap].at(j) == dcpbwt_build.phi->phi_inv_vec[hap].at(j));
    }
    for (unsigned int j = 0; j < dcpbwt_build.phi->phi_supp[hap].size(); ++j) {
      assert(dcpbwt_inserted.phi->phi_supp[hap].at(j) == dcpbwt_build.phi->phi_supp[hap].at(j));
    }
    for (unsigned int j = 0; j < dcpbwt_build.phi->phi_inv_supp[hap].size(); ++j) {
      assert(dcpbwt_inserted.phi->phi_inv_supp[hap].at(j) == dcpbwt_build.phi->phi_inv_supp[hap].at(j));
    }
    for (unsigned int j = 0; j < dcpbwt_build.phi->phi_supp_lcp[hap].size(); ++j) {
      assert(dcpbwt_inserted.phi->phi_supp_lcp[hap].at(j) == dcpbwt_build.phi->phi_supp_lcp[hap].at(j));
    }
  }
  std::cout << "PASSED!" << std::endl;
}

void Test_BottomUp_Delete(DCPBWT &dcpbwt) {
  cout << "Testing Deletion BottomUp i.e. from last haplotype to first haplotype in order...\n";
  int total = dcpbwt.M;
  for (int i = total - 1; i >= 0; --i) {
    dcpbwt.DeleteSingleHaplotype(dcpbwt.M - 1);
    cout << "Deleted hap " << i << "\n";
  }
  dcpbwt.DeleteSingleHaplotype(0);
}

void Test_TopDown_Delete(DCPBWT &dcpbwt) {
  cout << "Testing Deletion TopDown i.e. from first haplotype to last haplotype in order...\n";
  int total = dcpbwt.M;
  int hap_id = 0;
  for (int i = 0; i < total; ++i) {
    dcpbwt.DeleteSingleHaplotype(hap_id);
    cout << i << ". Deleted hap " << hap_id << "\n";
  }
  dcpbwt.DeleteSingleHaplotype(0);
}

void Insertion_Into_Empty_Panel(string &vcf_input, bool verbose) {
  DCPBWT dcpbwt;
  cout << "Inserting...\n";
  vector<vector<bool>> alleles;
  clock_t START_query_read = clock();
  ReadQueryVCF(vcf_input, alleles);
//  ReadQueryFile(vcf_input.c_str(), alleles);
  auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
  cout << "Time to read query alleles: " << time_read_query << " s\n";

  // go through all query haplotypes
  clock_t START_INSERT = clock();
  for (int i = 0; i < alleles.size(); ++i) {
    dcpbwt.InsertSingleHaplotype(alleles[i]);
    std::cout << "Inserted hap " << i << std::endl;
  }
  auto time_insert = (float) (clock() - START_INSERT) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";
  if (verbose){
    DCPBWT dcpbwt_build(vcf_input, verbose);
    Check_Insertion(dcpbwt, dcpbwt_build);
  }
}



void Insertion_Into_RefPanel(string &ref_vcf_input, string &query_vcf_input, bool verbose) {
  clock_t START = clock();
  DCPBWT dcpbwt(ref_vcf_input, verbose);
  auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
  cout << "Time to build: " << time_build << " secs.\n";

  cout << "Inserting...\n";
  vector<vector<bool>> alleles;
  clock_t START_query_read = clock();
  ReadQueryVCF(query_vcf_input, alleles);
//  ReadQueryFile(query_vcf_input.c_str(), alleles);
  auto time_read_query = (float) (clock() - START_query_read) / CLOCKS_PER_SEC;
  cout << "Time to read query alleles: " << time_read_query << " s\n";

  // go through all query haplotypes
  clock_t START_INSERT = clock();
  for(auto & allele : alleles) {
    dcpbwt.InsertSingleHaplotype(allele);
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
  std::string output_file;
  unsigned int length = 0;
  bool verbose = false;
  bool build = false;

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
      case 'b':build = true;
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
    DCPBWT dcpbwt(ref_vcf_input, verbose);
    dcpbwt.long_match_query(query_vcf_input, output_file, static_cast<unsigned int>(length), verbose, false);
    dcpbwt.PrintMemoryUsage(verbose);

  }else if (query){
    Insertion_Into_RefPanel(ref_vcf_input, query_vcf_input, verbose);
  } else {
    if (build){
      clock_t START = clock();
      DCPBWT dcpbwt(ref_vcf_input, verbose);
      auto time_build = (float)(clock() - START)/CLOCKS_PER_SEC;
      std::cout << "Build took " << time_build << " seconds." << std::endl;
      dcpbwt.PrintMemoryUsage(verbose);
    }else{
      Insertion_Into_Empty_Panel(ref_vcf_input, verbose);
    }
  }
  return 0;
}

