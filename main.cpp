#include <iostream>
#include <getopt.h>

#include "dcpbwt.h"
#include <random>

#include "utils.h"

void PrintHelp() {
  std::cout << "Usage: dcpbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input ref <path>\t vcf file for panel" << std::endl;
  std::cout << "  -q, --input query <path>\t vcf file for panel" << std::endl;
  std::cout << "  -l, --input length <int>\t of match" << std::endl;
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
  for (unsigned int i = 0; i < total; ++i) {
    cout << "Deleting hap " << i << "\n";
    dcpbwt.DeleteSingleHaplotype(0);
    cout << "Deleted hap " << i << "\n";
  }
  dcpbwt.DeleteSingleHaplotype(0);
}


void Test_Insertion(string &ref_vcf_input, string &query_vcf_input, bool verbose) {
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
  for (int i = 0; i < alleles.size(); ++i) {
    dcpbwt.InsertSingleHaplotype(alleles[i]);
//    cout << "Inserted query hap: " << i << "\n";
  }
  auto time_insert = (float) (clock() - START_INSERT) / CLOCKS_PER_SEC;
  cout << "Inserted " << alleles.size() << " haplotypes.\n";
  cout << "Insertion took: " << time_insert << " s.\n";

//  // go through all query haplotypes
//  clock_t START_DELETE = clock();
//  int total_delete = 1000;
//  for (int i = 0; i < total_delete; ++i) {
//    dcpbwt.DeleteSingleHaplotype(10);
//    cout << "Deleted query hap: " << i + 10 << "\n";
//  }
//  auto time_del = (float) (clock() - START_DELETE) / CLOCKS_PER_SEC;
//  cout << "Deleted " << total_delete << " haplotypes.\n";
//  cout << "Deletion took: " << time_del << " s.\n";
}

void Test_10Insertions(DCPBWT &dcpbwt) {
  vector<vector<bool>>
    data = {{false, true, false, false, true, false, true, false, false, false, true, true, true, false, true},
         {false, false, false, true, true, false, true, false, true, false, true, false, true, false, false},
         {true, true, false, false, true, false, true, true, false, false, false, true, false, true, true},
         {true, false, true, false, false, false, false, true, false, false, true, false, true, false, false},
         {false, false, true, false, true, false, true, true, true, false, true, true, false, true, true},
         {false, true, false, true, true, false, true, false, false, false, false, true, true, false, true},
         {true, true, false, false, false, false, true, false, false, false, true, true, false, true, false},
         {true, true, true, false, true, false, true, false, false, false, true, true, true, true, true},
         {true, true, false, false, false, false, true, false, true, false, true, false, false, false, false},
         {false, true, false, true, true, false, true, false, true, false, false, true, true, false, true}};

  for (int i = 0; i < data.size(); ++i) {
    dcpbwt.InsertSingleHaplotype(data[i]);
    cout << "Inserted Hap " << i << "\n";
  }

  assert(dcpbwt.M == 30);
  assert(dcpbwt.N == 15);
  vector<vector<int>> expected_pref_beg = {{0, 4, 17, 18, 22, 24, 26, 29},
                                           {4, 21, 25, 0, 17, 23, 26},
                                           {21, 24, 0, 23, 4, 18, 20, 27, 28},
                                           {21, 8, 9, 11, 16, 20, 25, 17},
                                           {8, 14, 20, 17, 22, 26, 24, 23, 18, 0, 1, 4, 25},
                                           {14, 17, 26, 4, 9},
                                           {14, 26, 23, 20, 18, 19, 1, 25, 17},
                                           {14, 23, 0, 1, 17, 22, 19},
                                           {14, 0, 9, 16, 12, 18, 28, 20, 21, 25, 29, 23, 24},
                                           {0, 8, 11},
                                           {0, 11, 18, 26, 19, 27, 25, 23, 1, 12, 29, 24},
                                           {0, 7, 25, 1, 2, 14, 29, 23, 12, 28, 24},
                                           {7, 1, 14, 9, 23, 28, 21, 18, 17, 4, 25, 2, 29, 11, 20, 12},
                                           {1, 28, 22, 11, 26, 12, 24, 8, 27},
                                           {28, 18, 11, 12, 23, 0, 26, 24},
                                           {28, 18}};

  vector<vector<int>> expected_pref_end = {{3, 16, 17, 21, 23, 25, 28, 29},
                                           {20, 24, 29, 3, 22, 23, 28},
                                           {21, 24, 3, 23, 16, 19, 26, 27, 28},
                                           {7, 8, 10, 15, 16, 20, 29, 27},
                                           {13, 15, 20, 17, 22, 28, 24, 23, 21, 0, 3, 16, 29},
                                           {15, 17, 0, 7, 29},
                                           {15, 28, 13, 24, 18, 21, 3, 29, 7},
                                           {15, 23, 18, 3, 20, 24, 29},
                                           {15, 0, 10, 11, 13, 26, 28, 27, 21, 25, 29, 22, 24},
                                           {16, 8, 24},
                                           {16, 11, 7, 20, 19, 27, 25, 23, 10, 21, 29, 8},
                                           {6, 19, 25, 1, 22, 10, 27, 23, 13, 21, 8},
                                           {19, 1, 15, 10, 23, 28, 16, 18, 17, 6, 25, 22, 29, 26, 27, 8},
                                           {10, 3, 22, 11, 26, 13, 24, 20, 27},
                                           {28, 3, 11, 15, 21, 22, 26, 27},
                                           {26, 27}};
  vector<vector<int>> expected_div_beg = {{0, 0, 0, 0, 0, 0, 0, 0},
                                           {1, 0, 0, 1, 0, 0, 0},
                                           {2, 0, 1, 0, 2, 0, 0, 0, 0},
                                           {3, 0, 0, 0, 0, 0, 0, 1},
                                           {4, 0, 0, 1, 0, 0, 3, 1, 2, 1, 0, 2, 0},
                                           {5, 1, 0, 2, 0},
                                           {6, 1, 3, 0, 2, 0, 1, 2, 6},
                                           {7, 3, 4, 4, 6, 1, 2},
                                           {8, 4, 2, 0, 0, 3, 0, 5, 4, 2, 0, 8, 3},
                                           {9, 5, 0},
                                           {10, 5, 3, 7, 3, 1, 4, 8, 5, 5, 2, 8},
                                           {11, 0, 4, 8, 0, 9, 7, 8, 9, 7, 8},
                                           {12, 8, 9, 4, 11, 9, 5, 5, 6, 4, 7, 8, 9, 11, 5, 9},
                                           {13, 11, 7, 11, 7, 9, 8, 10, 3},
                                           {14, 12, 11, 9, 11, 12, 11, 9},
                                           {15, 15}};
  cout << "############################## \n";
  cout << "Testing Prefix Samples Begin: \n";
  cout << "############################## \n";
  for (int col = 0; col <= dcpbwt.N; ++col) {
    assert(dcpbwt.columns[col].pref_samples_beg.size() == expected_pref_beg[col].size());
    for (int j = 0; j < expected_pref_beg[col].size(); ++j) {
      assert(dcpbwt.columns[col].pref_samples_beg.at(j) == expected_pref_beg[col][j]);
    }
    cout << "Col " << col << " passed !\n";
  }
  cout << "Passed !\n\n";

  cout << "############################## \n";
  cout << "Testing Prefix Samples End: \n";
  cout << "############################## \n";
  for (int col = 0; col <= dcpbwt.N; ++col) {
    assert(dcpbwt.columns[col].pref_samples_end.size() == expected_pref_end[col].size());
    for (int j = 0; j < expected_pref_end[col].size(); ++j) {
      assert(dcpbwt.columns[col].pref_samples_end.at(j) == expected_pref_end[col][j]);
    }
    cout << "Col " << col << " passed !\n";
  }
  cout << "Passed !\n\n";

  cout << "############################## \n";
  cout << "Testing Div Samples Beg: \n";
  cout << "############################## \n";
  for (int col = 0; col <= dcpbwt.N; ++col) {
    assert(dcpbwt.columns[col].div_samples_beg.size() == expected_div_beg[col].size());
    for (int j = 0; j < expected_div_beg[col].size(); ++j) {
      assert(dcpbwt.columns[col].div_samples_beg.at(j) == expected_div_beg[col][j]);
    }
    cout << "Col " << col << " passed !\n";
  }
  cout << "Passed !\n\n";

  cout << "############################## \n";
  cout << "Testing phi structures :\n";
  cout << "############################## \n";
  vector<vector<int>> expected_phi_supp = {
    {0, 29, 24, 21, 23, 15, 0, 0, 0, 21,   15}, //0
    {0, 21, 18, 23, 25, 19, 1,   20}, // 1
    {1, 25,   6}, // 2
    {2}, // 3
    {3, 4, 23, 3, 0, 17,     18}, // 4
    {4}, // 5
    {5}, // 6
    {6, 7, 8}, // 7
    {7, 8, 16, 24, 13}, // 8
    {8, 7, 0, 15, 1}, // 9
    {9 }, //10
    {10, 8, 16, 29, 22, 3, 28}, // 11
    {11, 10, 23, 27, 26, 11, 3}, // 12
    {12}, // 13
    {13, 14, 14, 14, 14, 22, 1, 19}, // 14
    {14}, // 15
    {15, 10, 0}, // 16
    {16, 3, 29, 20, 15, 29, 3, 18, 16}, // 17
    {17, 16, 23, 24, 13, 11, 16, 28, 26}, // 18
    {18, 24, 20, 7}, // 19
    {19, 16, 15, 13, 28, 26, 29}, // 20
    {20, 21, 21, 27, 28, 23}, // 21
    {21, 17, 20, 3, 10}, // 22
    {22, 3, 24, 28, 15, 29, 25, 27, 10, 15, 11}, // 23
    {23, 21, 28, 22, 29, 21, 13, 26, 22}, // 24
    {24, 20, 16, 3, 21, 27, 19, 6, 17}, // 25
    {25, 23, 22, 17, 15, 7, 11, 22, 21}, // 26
    {26, 19, 20, 24}, // 27
    {27, 26, 13, 23, 10, 28, 28}, // 28 last 28 might be excessive and we could may be add checks to not add that?
    {28, 25, 21, 10, 22, 25}, // 29
    };

  for(int i = 0; i < 30; ++i){
    assert(expected_phi_supp[i].size() == dcpbwt.phi->phi_supp[i].size());
    for(int j = 0; j < expected_phi_supp[i].size(); ++j){
      assert(expected_phi_supp[i][j] == dcpbwt.phi->phi_supp[i].at(j));
    }
    cout << "Passed phi supp for hap " << i << '\n';
  }
  cout << "\n";

  cout << "############################## \n";
  cout << "Testing phi inv structures :\n";
  cout << "############################## \n";
  vector<vector<int>> expected_phi_inv_supp = {
    {1, 4, 9, 16}, //0
    {2, 14, 9}, // 1
    {3}, // 2
    {4, 17, 23, 4, 25, 17, 22, 11, 12}, // 3
    {5}, // 4
    {6}, // 5
    {7, 25, 2}, // 6
    {8, 9, 7, 26, 19}, // 7
    {9, 11, 8, 8, 8, 7}, // 8
    {10}, // 9
    {11, 16, 12, 29, 23, 28, 22 }, //10
    {12, 18, 26, 12, 23}, // 11
    {13}, // 12
    {14, 20, 18, 28, 24, 8}, // 13
    {15}, // 14
    {16, 20, 17, 26, 23, 0, 9, 23, 0}, // 15
    {17, 18, 20, 25, 8, 11, 18, 17}, // 16
    {18, 22, 26, 4, 25}, // 17
    {19, 1, 17, 4}, // 18
    {20, 27, 25, 1, 14}, // 19
    {21, 25, 17, 22, 19, 27, 1}, // 20
    {22, 24, 0, 1, 25, 29, 24, 0, 26}, // 21
    {23, 26, 24, 14, 29, 11, 26, 24}, // 22
    {24, 26, 4, 18, 0, 1, 12, 28, 21}, // 23
    {25, 0, 23, 18, 19, 24, 24, 8, 27}, // 24
    {26, 29, 23, 1, 2, 29}, // 25
    {27, 28, 20, 12, 24, 18}, // 26
    {28, 27, 21, 25, 23, 12, 27, 27, 27}, // 27
    {29, 28, 28, 24, 23, 20, 21, 18, 11}, // 28 last 28 might be excessive and we could may be add checks to not add that?
    {29, 0, 17, 29, 29, 17, 29, 23, 24, 11, 20}, // 29
  };

  for(int i = 0; i < 30; ++i){
    assert(expected_phi_inv_supp[i].size() == dcpbwt.phi->phi_inv_supp[i].size());
    for(int j = 0; j < expected_phi_inv_supp[i].size(); ++j){
      assert(expected_phi_inv_supp[i][j] == dcpbwt.phi->phi_inv_supp[i].at(j));
    }
    cout << "Passed phi inv supp for hap " << i << '\n';
  }
}

void Test_Reverself(DCPBWT &dcpbwt) {
  unsigned int col = 0;
  unsigned int idx = 5;
  cout << "Col " << col << ", idx " << idx << " -> " << " Col " << col - 1 << ", idx " << dcpbwt.reverse_lf(col, idx)
       << "\n";
  col = 12;
  idx = 17;
  cout << "Col " << col << ", idx " << idx << " -> " << " Col " << col - 1 << ", idx " << dcpbwt.reverse_lf(col, idx)
       << "\n";
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
  if (!filesystem::exists(ref_vcf_input)){
    cerr << "Specified ref vcf : " << ref_vcf_input << " doesn't exist!\n";
    exit(EXIT_FAILURE);
  }
  if (!filesystem::exists(query_vcf_input)){
    cerr << "Specified query vcf : " << query_vcf_input << " doesn't exist!\n";
    exit(EXIT_FAILURE);
  }

//  Test_Insertion(ref_vcf_input, query_vcf_input, verbose);
//  Test_Deletion(ref_vcf_input, query_vcf_input, verbose);
  DCPBWT dcpbwt(ref_vcf_input, verbose);
  dcpbwt.DeleteSingleHaplotype(19);
  dcpbwt.DeleteSingleHaplotype(4);
  vector<bool> query = {false, true, true, false, true, false, true, false, false, false, false, false, true, false, true};
  vector<pair<unsigned int, unsigned int>> insertion_indices(dcpbwt.N + 1); // stores {index, hapid}
  insertion_indices[0].first = 4;
  insertion_indices[0].second = 5;

  // Calculate insertion indices
  for (unsigned int col = 0; col < query.size(); ++col) {
    insertion_indices[col + 1] = dcpbwt.w_mod(insertion_indices[col].first, col, query[col], insertion_indices[col].second);
  }

  for(int i = 0; i < insertion_indices.size(); ++i){
    cout << insertion_indices[i].first << " : " << insertion_indices[i].second << "\n";
  }
//  Test_10Insertions(dcpbwt);
  return (EXIT_SUCCESS);
}
