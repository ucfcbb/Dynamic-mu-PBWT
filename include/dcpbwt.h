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
    std::vector<unsigned int> num_zeros;


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

    [[nodiscard]] unsigned int get_run_idx(const unsigned int col, const unsigned int i) const {
      return combined[col].search(i+1);
    }

    // returns the head of a run
    [[nodiscard]] unsigned int get_run_head(const unsigned int col, const unsigned int run_idx) const {
      if (run_idx >= combined[col].size()) {
        cerr << "Out of bounds: Accessing run that doesn't exist!\n";
        exit(EXIT_FAILURE);
      }
      if (run_idx == 0) {
        return 0;
      }
      return combined[col].psum(run_idx - 1);
    }

    bool get_run_val(const unsigned int col, const unsigned int run_idx) {
      if (run_idx % 2 == 0) {
        if (start_with_zero[col]) {
          return false;
        }
        return true;
      }
      if (start_with_zero[col]) {
        return true;
      }
      return false;
    }

    unsigned int lf(const unsigned int col, const unsigned int i, const bool value) {
      if (i == M) {
        if (value) // val is 1
          return M;
        return num_zeros[col]; // val is 0
      }

      unsigned int run_idx = get_run_idx(col, i);
      unsigned int run_head = get_run_head(col, run_idx);
      unsigned int offset = i - run_head;

      if (value != get_run_val(col, run_idx)) {
        offset = 0;
      }

      auto uv  = uv_trick(col, i);

      if (value) {
        return num_zeros[col] + uv.second + offset;
      }
      return uv.first + offset;
    }


    unsigned int zeros_before(const unsigned int col, const unsigned int run_idx) {
      unsigned int retval =0;
      assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
      if (start_with_zero[col]) {
        if (run_idx % 2 == 0)
         retval = zeros[col].psum((run_idx - 2)/2);
        else
         retval = zeros[col].psum((run_idx - 1)/2);
      } else {
        if (run_idx % 2 == 0)
         retval = zeros[col].psum((run_idx - 1)/2);
        else
         retval = zeros[col].psum((run_idx - 2)/2);
      }
      return retval;
    }

    unsigned int ones_before(const unsigned int col, const unsigned int run_idx) {
      unsigned int retval =0;
      assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
      if (start_with_zero[col]) {
        if (run_idx % 2 == 0)
         retval = ones[col].psum((run_idx - 1)/2);
        else
         retval = ones[col].psum((run_idx - 2)/2);
      } else {
        if (run_idx % 2 == 0)
         retval = ones[col].psum((run_idx - 2)/2);
        else
         retval = ones[col].psum((run_idx - 1)/2);
      }
      return retval;
    }

    // returns the number of zeros and ones before ith haplotype in a column
    pair<unsigned int, unsigned int> uv_trick(const unsigned int col, const unsigned i) {
      unsigned int run_idx = get_run_idx(col, i);
      unsigned int u = 0, v = 0;
      if (run_idx == 0) {
        return make_pair(0, 0);
      }
      if (run_idx == 1) {
        if(start_with_zero[col]) {
          return make_pair(zeros_before(col, run_idx), 0);
        }
        return make_pair(0, ones_before(col, run_idx));
      }
      u = zeros_before(col, run_idx);
      v = ones_before(col, run_idx);
      return make_pair(u, v);
    }

    void InsertSinglelHaplotype(std::vector<bool>& query) {
      assert(query.size() == N);
      unsigned int idx = M;
      cout << "At col0: idx : " << idx << "\n";
      for(unsigned int col = 1; col < query.size(); ++col) {
        idx = lf(col-1, idx, query[col-1]);
        cout << "At col" << col << " idx: " << idx << "\n";
        // Insert(col, idx, query[col]);
        // update all auxiliary data structures
      }
    }

    // Build the reference panel
    void Build (std::vector<std::vector<bool>>& alleles){
      // TODO: NEED TO implement UPDATE Phi data structures

      vector<int> u, v;
      vector<int> freq;
      int total_runs = 0;
      int col = 0;
      int cnt = 1;
      bool prev_allele = false;
      vector<int> prefix_arr(M, 0);
      std::iota(prefix_arr.begin(), prefix_arr.end(), 0);

      while (col < N) {
        packed_spsi temp_zeros;
        packed_spsi temp_ones;
        packed_spsi temp_combined;
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
        num_zeros.push_back(u.size());

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
