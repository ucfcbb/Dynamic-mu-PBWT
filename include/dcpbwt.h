#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <sstream>

#include "dynamic/dynamic.hpp"
#include "phi.h"
#include "dcpbwt_column.h"

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
    std::vector<dcpbwt_column> columns;
    phi_ds* phi=nullptr;

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
      return columns[col].combined.search(i+1);
    }

    // returns the head of a run
    [[nodiscard]] unsigned int get_run_head(const unsigned int col, const unsigned int run_idx) const {
      if (run_idx >= columns[col].combined.size()) {
        cerr << "Out of bounds: Accessing run that doesn't exist!\n";
        exit(EXIT_FAILURE);
      }
      if (run_idx == 0) {
        return 0;
      }
      return columns[col].combined.psum(run_idx - 1);
    }

    [[nodiscard]] bool get_run_val(const unsigned int col, const unsigned int run_idx) const {
      if (run_idx % 2 == 0) {
        if (columns[col].start_with_zero) {
          return false;
        }
        return true;
      }
      if (columns[col].start_with_zero) {
        return true;
      }
      return false;
    }

    unsigned int lf(const unsigned int col, const unsigned int i, const bool value) {
      if (i == M) {
        if (value) // val is 1
          return M;
        return columns[col].num_zeros; // val is 0
      }

      unsigned int run_idx = get_run_idx(col, i);
      unsigned int run_head = get_run_head(col, run_idx);
      unsigned int offset = i - run_head;

      if (value != get_run_val(col, run_idx)) {
        offset = 0;
      }

      auto uv  = uv_trick(col, i);

      if (value) {
        return columns[col].num_zeros + uv.second + offset;
      }
      return uv.first + offset;
    }

    unsigned int zeros_before(const unsigned int col, const unsigned int run_idx) {
      unsigned int retval =0;
      assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
      if (columns[col].start_with_zero) {
        if (run_idx % 2 == 0)
         retval = columns[col].zeros.psum((run_idx - 2)/2);
        else
         retval = columns[col].zeros.psum((run_idx - 1)/2);
      } else {
        if (run_idx % 2 == 0)
         retval = columns[col].zeros.psum((run_idx - 1)/2);
        else
         retval = columns[col].zeros.psum((run_idx - 2)/2);
      }
      return retval;
    }

    unsigned int ones_before(const unsigned int col, const unsigned int run_idx) {
      unsigned int retval =0;
      assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
      if (columns[col].start_with_zero) {
        if (run_idx % 2 == 0)
         retval = columns[col].ones.psum((run_idx - 1)/2);
        else
         retval = columns[col].ones.psum((run_idx - 2)/2);
      } else {
        if (run_idx % 2 == 0)
         retval = columns[col].ones.psum((run_idx - 2)/2);
        else
         retval = columns[col].ones.psum((run_idx - 1)/2);
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
        if(columns[col].start_with_zero) {
          return make_pair(zeros_before(col, run_idx), 0);
        }
        return make_pair(0, ones_before(col, run_idx));
      }
      u = zeros_before(col, run_idx);
      v = ones_before(col, run_idx);
      return make_pair(u, v);
    }

    void Insert(const unsigned int col, const unsigned i, const bool allele) {
      // Inserting allele at the bottom of a column
      if (i == M) {
        if (get_run_val(col, get_run_idx(col, i-1)) != allele) {
          columns[col].combined.push_back(1);
          if (allele) {
            columns[col].ones.push_back(1);
          } else {
            columns[col].zeros.push_back(1);
            ++columns[col].num_zeros;
          }
          // insert new run head at the bottom
          columns[col].pref_samples_beg.push_back(M);
          columns[col].pref_samples_end.push_back(M);

          // update phi structure
          this->phi->phi_vec[M].set(col, true);
          this->phi->phi_inv_vec[M].set(col, true);
        } else {
          columns[col].combined.increment(columns[col].combined.size()-1, 1);
          if (allele) {
            columns[col].ones.increment(columns[col].ones.size() - 1, 1);
          } else {
            columns[col].zeros.increment(columns[col].zeros.size() - 1, 1);
            ++columns[col].num_zeros;
          }
          // unset the bit for previous haplotype that was at the end of run in this column
          this->phi->phi_inv_vec[columns[col].pref_samples_end[columns[col].pref_samples_end.size()-1]].set(col, false);
          // replace the pref val at end of a run for column 'col' to M
          columns[col].pref_samples_end.set(columns[col].pref_samples_end.size() - 1, M);
          // update phi structure
          this->phi->phi_inv_vec[M].set(col, true);
        }
        return;
      }

      // Inserting allele at the top of a column
      if (i == 0) {
        if (columns[col].start_with_zero) {
          if (allele) { // insert 1
            columns[col].combined.insert(i, 1);
            columns[col].ones.insert(i, 1);
            columns[col].start_with_zero = false;

            columns[col].pref_samples_beg.insert(i, M);
            columns[col].pref_samples_end.insert(i, M);

            // phi support
            this->phi->phi_vec[M].set(col, true);
            this->phi->phi_inv_vec[M].set(col, true);

          } else { // insert 0
            columns[col].combined.increment(i, 1);
            columns[col].zeros.increment(i, 1);
            ++columns[col].num_zeros;

            // update(unset) old sample beg
            auto curr_top_hap = columns[col].pref_samples_beg.at(0);
            this->phi->phi_vec[curr_top_hap].set(col, false);

            // new sample beg
            columns[col].pref_samples_beg.set(0, M);
            this->phi->phi_vec[M].set(col, true);
          }
        } else {
          if (allele) {
            columns[col].combined.increment(i, 1);
            columns[col].ones.increment(i, 1);

            // update(unset) old sample beg
            auto old_top_hap = columns[col].pref_samples_beg.at(0);
            this->phi->phi_vec[old_top_hap].set(col, false);

            // new sample beg
            columns[col].pref_samples_beg.set(0, M);
            this->phi->phi_vec[M].set(col, true);

          } else {
            columns[col].combined.insert(i, 1);
            columns[col].zeros.insert(i, 1);
            columns[col].start_with_zero = true;
            ++columns[col].num_zeros;

            columns[col].pref_samples_beg.insert(i, M);
            columns[col].pref_samples_end.insert(i, M);

            // phi support
            this->phi->phi_vec[M].set(col, true);
            this->phi->phi_inv_vec[M].set(col, true);
          }
        }
        return;
      }

      // get run index and run value
      unsigned int run_idx = get_run_idx(col, i);
      bool run_value = get_run_val(col, run_idx);
      unsigned int run_head = get_run_head(col, run_idx);

      // for insertions at the boundary of runs
      if (run_head == i) {
        // update previous run's end boundaries
        if (run_value != allele) {
          --run_idx;
          columns[col].combined.increment(run_idx, 1);
          if(allele) {
            columns[col].ones.increment(run_idx/2, 1);
          } else {
            columns[col].zeros.increment(run_idx/2, 1);
            ++columns[col].num_zeros;
          }
          // since allele is being added at tne end of prev. run

          //unset old hap
          auto old_end_hap = columns[col].pref_samples_end.at(run_idx);
          this->phi->phi_inv_vec[old_end_hap].set(col, false);

          // set the new one
          columns[col].pref_samples_end.set(run_idx, M);
          this->phi->phi_inv_vec[M].set(col, true);
        } else {
          // change to new run head
          columns[col].combined.increment(run_idx, 1);
          if(allele) {
            columns[col].ones.increment(run_idx/2, 1);
          } else {
            columns[col].zeros.increment(run_idx/2, 1);
            ++columns[col].num_zeros;
          }
          // update(unset) old sample beg
          auto old_top_hap = columns[col].pref_samples_beg.at(run_idx);
          this->phi->phi_vec[old_top_hap].set(col, false);

          // set the new sample beg
          columns[col].pref_samples_beg.set(run_idx, M);
          this->phi->phi_vec[M].set(col, true);
        }
        return;
      }

      // if run value is same as the allele
      if (run_value == allele) {
        columns[col].combined.increment(run_idx, 1);
        if(allele) {
          columns[col].ones.increment(run_idx/2, 1);
        } else {
          columns[col].zeros.increment(run_idx/2, 1);
          ++columns[col].num_zeros;
        }
        return;
      }


      // split insertion in the middle
      unsigned int left_rem = i - run_head;
      unsigned int right_rem = columns[col].combined.at(run_idx) - left_rem;
      columns[col].combined.set(run_idx, left_rem);
      columns[col].combined.insert(run_idx + 1, 1);
      columns[col].combined.insert(run_idx + 2, right_rem);
      unsigned int  new_idx = run_idx/2 + 1;

      // TODO
      // How to find the sample beg of the right half of the split run i.e.
      // Alleles: 1 1 0 0 0 (1) 0 0 0 1 1 1
      // PrefBeg: 3   5      x  ?     8
      // perhaps need to implement phi first and then use phi to get that pref value


      // find hapID that is before and after inserted haplotype
      int loop_cnt = i - run_head - 1;
      unsigned int hap_before = 0;
      unsigned int start_hap = columns[col].pref_samples_beg.at(run_idx);
      for(int j = 0; j < loop_cnt;++j) {
        hap_before = this->phi->phi_inv(start_hap, col).value();
        start_hap = hap_before;
      }
      unsigned int hap_after = this->phi->phi_inv(hap_before, col).value();

      // insert left to right
      columns[col].pref_samples_beg.insert(run_idx + 1, M);
      columns[col].pref_samples_beg.insert(run_idx + 2, hap_after);

      // insert right to left
      columns[col].pref_samples_end.insert(run_idx, M);
      columns[col].pref_samples_end.insert(run_idx, hap_before);

      // update phi structures
      // Inserted Haplotype
      this->phi->phi_vec[M].set(col, true);
      this->phi->phi_inv_vec[M].set(col, true);

      // Haplotypes that are before and after the inserted haplotype
      this->phi->phi_vec[hap_after].set(col, true);
      this->phi->phi_inv_vec[hap_before].set(col, true);

      if (allele) {
        assert(new_idx <= columns[col].ones.size());
        columns[col].ones.insert(new_idx, 1);
      } else {
        assert(new_idx <= columns[col].zeros.size());
        columns[col].zeros.insert(new_idx, 1);
        ++columns[col].num_zeros;
      }
    }

    void InsertSinglelHaplotype(std::vector<bool>& query) {
      assert(query.size() == N);
      vector<unsigned int> insertion_indices(N,0);
      insertion_indices[0] = M;

      // calculate insertion indices
      for(unsigned int col = 1; col < query.size(); ++col) {
        insertion_indices[col] = lf(col-1, insertion_indices[col-1], query[col-1]);
      }
      // initialize new phi structure
      suc_bv tmp_b, tmp_e;
      for(unsigned int col = 0; col < query.size(); ++col) {
        tmp_b.push_back(false);
        tmp_e.push_back(false);
      }
      this->phi->phi_vec.push_back(tmp_b);
      this->phi->phi_inv_vec.push_back(tmp_e);
      // TODO: how to update the support vectors?

      // perform insertion
      for(unsigned int col = 0; col < query.size(); ++col) {
        Insert(col, insertion_indices[col], query[col]);
      }

      // update support vectors
      // for (unsigned int col = 0; col < this->N; col++) {
      //   for (unsigned int j = 0; j < columns[col].pref_samples_beg.size(); j++) {
      //     // support phi panel (if we are in the first run we use default
      //     // value)
      //     if (j == 0) {
      //       this->phi->phi_supp[columns[col].pref_samples_beg[j]].push_back(UINT_MAX);
      //     } else {
      //       this->phi->phi_supp[columns[col].pref_samples_beg[j]].push_back(columns[col].pref_samples_end[j - 1]);
      //     }
      //     // use sample end to compute phi_inv panel
      //     if (j == columns[col].pref_samples_beg.size() - 1) {
      //       this->phi->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
      //               UINT_MAX);
      //     } else {
      //       this->phi->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
      //               columns[col].pref_samples_beg[j+1]);
      //     }
      //   }
      // }

      // increment total # of haplotypes
      ++M;
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
      // TODO: Possible optimization using sdsl int vec
      vector<int> prefix_arr(M, 0);
      std::iota(prefix_arr.begin(), prefix_arr.end(), 0);

      // TODO: Possible optimization using sdsl int vec
      vector<vector<unsigned int>> sites_where_sample_beg(M);
      vector<vector<unsigned int>> sites_where_sample_end(M);

      while (col < N) {
        packed_spsi temp_zeros;
        packed_spsi temp_ones;
        packed_spsi temp_combined;
        packed_spsi temp_sample_beg;
        packed_spsi temp_sample_end;
        bool start_with_zero = false;

        for (int i = 0; i < M; ++i) {
          // first allele
          if (i == 0) {
            if (alleles[col][prefix_arr[i]]) { // allele: 1
              v.push_back(prefix_arr[i]);
            } else { // allele: 0
              u.push_back(prefix_arr[i]);
              start_with_zero = true;
            }
            temp_sample_beg.push_back(prefix_arr[i]);
            prev_allele = alleles[col][prefix_arr[i]];
            cnt = 1;
            continue;
          }

          if (alleles[col][prefix_arr[i]] != prev_allele) {
            freq.push_back(cnt);
            prev_allele = alleles[col][prefix_arr[i]];
            cnt = 1;
            temp_sample_beg.push_back(prefix_arr[i]);
            temp_sample_end.push_back(prefix_arr[i - 1]);
          } else {
            ++cnt;
          }

          if (alleles[col][prefix_arr[i]]) {
            v.push_back(prefix_arr[i]);
          } else {
            u.push_back(prefix_arr[i]);
          }
        }
        temp_sample_end.push_back(prefix_arr[M-1]);

        // edge case
        if (cnt > 0)
          freq.push_back(cnt);

        // populate dynamic data structures
        for(int i = 0; i < freq.size(); ++i) {
          temp_combined.push_back(freq[i]);
          if (start_with_zero) {
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

        assert(temp_combined.size() == temp_sample_beg.size());
        assert(temp_combined.size() == temp_sample_end.size());

        for(auto i = 0; i < temp_sample_beg.size(); ++i) {
          sites_where_sample_beg[temp_sample_beg.at(i)].push_back(col);
          sites_where_sample_end[temp_sample_end.at(i)].push_back(col);
        }

        // build each column
        dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                           (std::move(temp_sample_beg)), (std::move(temp_sample_end)),
                           start_with_zero, u.size());
        columns.emplace_back(coln);


        total_runs += freq.size();
        // next col prefix arr
        prefix_arr.clear();
        prefix_arr.insert(prefix_arr.end(), u.begin(), u.end());
        prefix_arr.insert(prefix_arr.end(), v.begin(), v.end());
        ++col;
        u.clear();
        v.clear();
        freq.clear();
      } // run-through all sites

      assert(columns.size() == N);
      // build phi data-structure
      this->phi = new phi_ds(columns, M, N, sites_where_sample_beg, sites_where_sample_end, prefix_arr, false);
      assert(col == N);
      total_runs += freq.size();
      cout << "Phi support size (in bits) = " << this->phi->size_in_bytes(false) << "\n";
      cout << "Avg runs = " << static_cast<float>(total_runs)/N << "\n";
    }
};
