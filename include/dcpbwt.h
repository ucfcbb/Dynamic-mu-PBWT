#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <sstream>

#include "dynamic/dynamic.hpp"

using namespace dyn;
using namespace std;

//typedef spsi<packed_vector, 256, 2> my_spsi;

void ReadVCF(string& filename, vector<vector<bool>>& alleles){
  int M = 0;
  int N = 0;
  std::string line= "##";
  std::ifstream inFile(filename);
  if (inFile.is_open()){
    while(line[1] == '#'){
      getline(inFile, line);
    }

    // read individual ids
    std::istringstream iss(line);
    std::string token;
    for(int i = 0; i <9; ++i){
      iss >> token;
    }
    while(iss >> token){
      M += 2;
      // individual_ids.push_back(token);
      continue;
    }

    // go through all sites
    while(getline(inFile, line)){
      ++N;
      std::istringstream iss(line);
      token = "";
      std::vector<bool> single_col;
      for(int i = 0; i < (M/2) +9 ; ++i){
        iss >> token;
        if (i < 9){
          continue;
        }
        single_col.push_back(static_cast<bool>(token[0] - '0'));
        single_col.push_back(static_cast<bool>(token[2] - '0'));
      }
      alleles.push_back(single_col);
    }
    inFile.close();
  }else{
    std::cerr << "Couldn't find : " << filename << "\n";
    exit(1);
  }
  std::cout << "M (# of haplotypes) = " << M << " : N (# of sites) = " << N << "\n";
}

class DCPBWT{
  public:
    unsigned int M; // #haplotypes
    unsigned int N; // #sites
    std::vector<packed_spsi> zeros;
    std::vector<packed_spsi> ones;
    std::vector<packed_spsi> combined;
    std::vector<bool> start_with_zero;



    // constructor
    DCPBWT(std::string ref_vcf_file, bool verbose){
      std::vector<std::vector<bool>> alleles;
      // extract alleles from VCF
      ReadVCF(ref_vcf_file, alleles);
      N = alleles.size();
      M = alleles[0].size();
      Build(alleles);

      // build the ref panel
    }

    [[nodiscard]] unsigned int run_idx(const unsigned int col, const unsigned int i) const {
      return combined[col].search(i+1);
    }

    // returns the head of a run
    unsigned int run_head(const unsigned int col, const unsigned int run_idx) {
      if (run_idx >= combined[col].size()) {
        cerr << "Out of bounds: Accessing run that doesn't exist!\n";
        exit(EXIT_FAILURE);
      }
      if (run_idx == 0) {
        return 0;
      }
      return combined[col].psum(run_idx - 1);
    }

    // returns the number of zeros and ones before ith haplotype in a column
    pair<unsigned int, unsigned int> uv_trick(const unsigned int col, const unsigned i) {

    }




    // Build the reference panel
    void Build (std::vector<std::vector<bool>>& alleles){
      // TODO: NEED TO implement UPDATE Phi data structures
      packed_spsi temp_zeros;
      packed_spsi temp_ones;
      packed_spsi temp_combined;

      vector<int> u, v;
      vector<int> freq;
      int total_runs = 0;
      int col = 0;
      int cnt = 1;
      bool prev_allele = false;
      vector<int> prefix_arr(M, 0);
      std::iota(prefix_arr.begin(), prefix_arr.end(), 0);

      while (col < N) {
        for (int i = 0; i < M; ++i) {
          if (i == 0) {
            if (alleles[col][prefix_arr[i]]) { // allele: 1
              v.push_back(prefix_arr[i]);
              start_with_zero.push_back(false);
            } else { // allele: 0
              u.push_back(prefix_arr[i]);
              start_with_zero.push_back(true);
            }
            prev_allele = alleles[col][prefix_arr[i]];
            cnt = 1;
            continue;
          }

          if (alleles[col][prefix_arr[i]] != prev_allele) {
            freq.push_back(cnt);
            prev_allele = alleles[col][prefix_arr[i]];
            cnt = 1;
          } else {
            ++cnt;
          }
          if (alleles[col][prefix_arr[i]]) {
            v.push_back(prefix_arr[i]);
          } else {
            u.push_back(prefix_arr[i]);
          }
        }
        // edge case
        if (cnt > 0)
          freq.push_back(cnt);

        // populate dynamic data structures
        for(int i = 0; i < freq.size(); ++i) {
          temp_combined.push_back(freq[i]);
          if (start_with_zero[col]) {
            if (i % 2 == 0) {
              temp_zeros.push_back(freq[i]);
            } else {
              temp_ones.push_back(freq[i]);
            }
          } else {
            if (i % 2 == 0) {
              temp_ones.push_back(freq[i]);
            } else {
              temp_zeros.push_back(freq[i]);
            }
          }
        }

        combined.push_back(temp_combined);
        zeros.push_back(temp_zeros);
        ones.push_back(temp_ones);

        total_runs += freq.size();
        // next col prefix arr
        prefix_arr.clear();
        prefix_arr.insert(prefix_arr.end(), u.begin(), u.end());
        prefix_arr.insert(prefix_arr.end(), v.begin(), v.end());
        ++col;
        u.clear();
        v.clear();
        freq.clear();
      }
      assert(col == N);
      total_runs += freq.size();
      cout << "Avg runs = " << static_cast<float>(total_runs)/N << "\n";
    }
};
