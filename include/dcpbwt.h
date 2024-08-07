#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "dynamic/dynamic.hpp"

using namespace dyn;

//typedef spsi<packed_vector, 256, 2> my_spsi;

class DCPBWT{
  public:
    const unsigned int M; // #haplotypes
    const unsigned int N; // #sites
    std::vector<packed_spsi> zeros;
    std::vector<packed_spsi> ones;
    std::vector<packed_spsi> combined;
    std::vector<bool> start_with_zero;

    // constructor
    DCPBWT(std::string ref_vcf_file, bool verbose){
      std::vector<std::vector<bool>> alleles;
      // extract alleles from VCF
      ReadFromVCF(ref_vcf_file, alleles);
      N = alleles.size();
      M = alleles[0].size();

      // build the ref panel
    }

    Build (std::vector<std::vector<bool>>& alleles){
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
            heads.push_back(i);
            continue;
          }

          if (alleles[col][prefix_arr[i]] != prev_allele) {
            if (alleles[col][prefix_arr[i]]){
              temp_zeros.push_back(u.size());
            } else {
              temp_ones.push_back(v.size());
            }

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

        //cout << "Col " << col << ", #runs = " << freq.size() << "\n\n";

        // insert into B+-tree
        columns[col] = new my_spsi();
        for(const int fq: freq) {
          columns[col].push_back(fq);
        }

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
      cout << "Avg runs = " << (float)total_runs/N << "\n";
    }
};
