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

inline pair<unsigned int, unsigned int> ReadVCF(string& filename) {
  unsigned int M = 0;
  unsigned int N = 0;
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
    }
    inFile.close();
  }else{
    std::cerr << "Couldn't find : " << filename << "\n";
    exit(1);
  }
  std::cout << "M (# of haplotypes) = " << M << " : N (# of sites) = " << N << "\n";
  return {M, N};
}

inline void ReadVCF(string& filename, vector<vector<bool>>& alleles){
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
      // std::vector<std::vector<bool>> alleles;
      // extract alleles from VCF
      auto retval = ReadVCF(ref_vcf_file);
      this->M = retval.first;
      this->N = retval.second;

      // build the ref panel
      // Build(alleles);
      BuildFromVCF(ref_vcf_file, verbose);
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

    std::pair<unsigned int,unsigned int> column_end_lf(unsigned int col, bool val) {
      if (val) { return {this->M, UINT_MAX};} // last row, invalid hapID
      unsigned int ind = columns[col].num_zeros;
      unsigned int ridx = 0; // assume column starts with a run of ones
      if (columns[col].start_with_zero){ // the second run at index 1 will be run of ones
        ridx = 1;
      }
      return {ind, columns[col].pref_samples_beg.at(ridx)};
    }

     std::pair<int,int> w_mod(const unsigned int i, const unsigned int col, const bool val, const unsigned int pref_val){
        if (i == this->M) {return column_end_lf(col, val);}
        const unsigned int run_idx = get_run_idx(col, i);

        if (const bool rval = get_run_val(col, run_idx); rval == val){
          return {lf(col, i, val), pref_val};
        }

        if (run_idx + 1 >= columns[col].combined.size())
          return column_end_lf(col, val);
        return {lf(col, get_run_head(col, run_idx + 1), val), columns[col].pref_samples_beg[run_idx+1]}; // sample_beg fetch pref-val
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

    void Insert(const unsigned int col, const unsigned i, const unsigned int hap_id, bool allele, packed_spsi& temp_supp, packed_spsi& temp_inv_supp) {
      // assumes Non-empty ref panel (i.e. a panel is built prior to insertion)
      // Inserting allele at the bottom of a column
      if (i == M) {
        // inserting a new run
        if (get_run_val(col, get_run_idx(col, i-1)) != allele) {
          columns[col].combined.push_back(1);
          if (allele) {
            columns[col].ones.push_back(1);
          } else {
            columns[col].zeros.push_back(1);
            ++columns[col].num_zeros;
          }

          // get old "end of run" hap and update it's bottom
          assert (columns[col].pref_samples_end.size() > 0);
          auto old_bottom_hap = columns[col].pref_samples_end.at(columns[col].pref_samples_end.size() - 1);
          auto col_rank = this->phi->phi_inv_vec[old_bottom_hap].rank1(col);
          this->phi->phi_inv_supp[old_bottom_hap].set(col_rank, M);

          // insert new run head at the bottom
          columns[col].pref_samples_beg.push_back(M);
          columns[col].pref_samples_end.push_back(M);

          // update phi structure
          this->phi->phi_vec[M].set(col, true);
          this->phi->phi_inv_vec[M].set(col, true);

          temp_supp.push_back(old_bottom_hap);
          std::cout << "size of temp-supp = " << temp_supp.size() << " , at 0 = " << temp_supp.at(0) << "\n";
          // UINT_MAX represents nothing below it
          temp_inv_supp.push_back(UINT_MAX);
        } else { // modifying end of existing run
          columns[col].combined.increment(columns[col].combined.size()-1, 1);
          if (allele) {
            columns[col].ones.increment(columns[col].ones.size() - 1, 1);
          } else {
            columns[col].zeros.increment(columns[col].zeros.size() - 1, 1);
            ++columns[col].num_zeros;
          }
          // unset the bit for previous haplotype that was at the end of run in this column
          auto old_bottom_hap = columns[col].pref_samples_end.at(columns[col].pref_samples_end.size() - 1);
          auto col_rank = this->phi->phi_inv_vec[old_bottom_hap].rank1(col);
          this->phi->phi_inv_supp[old_bottom_hap].remove(col_rank);
          this->phi->phi_inv_vec[old_bottom_hap].set(col, false);

          // replace the pref val at end of a run for column 'col' to M
          columns[col].pref_samples_end.set(columns[col].pref_samples_end.size() - 1, M);
          // update phi structure
          this->phi->phi_inv_vec[M].set(col, true);
          temp_inv_supp.push_back(UINT_MAX);
        }
        return;
      }

      // Inserting allele at the top of a column
      if (i == 0) {
        if (columns[col].start_with_zero) {
          if (allele) { // insert 1
            // TODO: refactor to insert_new_run();
            columns[col].combined.insert(i, 1);
            columns[col].ones.insert(i, 1);
            columns[col].start_with_zero = false;

            // update for the haplotype that's at the top
            auto old_top_hap = hap_id;
            this->phi->phi_supp[old_top_hap].set(i, M);

            // update the data-structures for the inserted haplotype
            // update run beg/end pref samples
            columns[col].pref_samples_beg.insert(i, M);
            columns[col].pref_samples_end.insert(i, M);

            // update phi support data structures
            this->phi->phi_vec[M].set(col, true);
            this->phi->phi_inv_vec[M].set(col, true);
            // nothing above inserted hap
            temp_supp.push_back(UINT_MAX);
            temp_inv_supp.push_back(old_top_hap);

          } else { // insert 0
            columns[col].combined.increment(i, 1);
            columns[col].zeros.increment(i, 1);
            ++columns[col].num_zeros;

            // unset old sample beg
            auto old_top_hap = hap_id;
            auto col_rank = this->phi->phi_vec[old_top_hap].rank1(col);
            this->phi->phi_supp[old_top_hap].remove(col_rank);
            this->phi->phi_vec[old_top_hap].set(col, false);

            // update new sample beg
            columns[col].pref_samples_beg.set(i, M);
            this->phi->phi_vec[M].set(col, true);
            temp_supp.push_back(UINT_MAX);
          }
        } else {
          if (allele) {
            columns[col].combined.increment(i, 1);
            columns[col].ones.increment(i, 1);

            // remove & unset old sample beg
            auto old_top_hap = hap_id;
            auto col_rank = this->phi->phi_vec[old_top_hap].rank1(col);
            this->phi->phi_supp[old_top_hap].remove(col_rank);
            this->phi->phi_vec[old_top_hap].set(col, false);

            // insert new sample beg
            columns[col].pref_samples_beg.set(i, M);
            this->phi->phi_vec[M].set(col, true);
            temp_supp.push_back(UINT_MAX);
          } else {
            // insert new run
            columns[col].combined.insert(i, 1);
            columns[col].zeros.insert(i, 1);
            columns[col].start_with_zero = true;
            ++columns[col].num_zeros;

            // update old top's related data-structures
            auto old_top_hap = hap_id;
            auto col_rank = this->phi->phi_vec[old_top_hap].rank1(col);
            this->phi->phi_supp[old_top_hap].set(col_rank , M);

            columns[col].pref_samples_beg.insert(i, M);
            columns[col].pref_samples_end.insert(i, M);

            // phi support
            this->phi->phi_vec[M].set(col, true);
            this->phi->phi_inv_vec[M].set(col, true);
            temp_supp.push_back(UINT_MAX);
            temp_inv_supp.push_back(old_top_hap);
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
          // since allele is being added at the end of prev. run
          // unset old hap
          auto old_end_hap = columns[col].pref_samples_end.at(run_idx);
          auto col_rank = this->phi->phi_inv_vec[old_end_hap].rank1(col);
          this->phi->phi_inv_supp[old_end_hap].remove(col_rank);
          this->phi->phi_inv_vec[old_end_hap].set(col, false);

          // update the run_idx + 1's top as well
          col_rank = this->phi->phi_vec[hap_id].rank1(col);
          this->phi->phi_supp[hap_id].set(col_rank, M);

          // set the new one
          columns[col].pref_samples_end.set(run_idx, M);
          this->phi->phi_inv_vec[M].set(col, true);
          temp_inv_supp.push_back(hap_id);
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
          auto old_top_hap = hap_id;
          auto col_rank = this->phi->phi_vec[old_top_hap].rank1(col);
          this->phi->phi_supp[old_top_hap].remove(col_rank);
          this->phi->phi_vec[old_top_hap].set(col, false);

          // update what will be below the previous run's end
          assert(run_idx - 1 < UINT_MAX);
          auto above_old_top_hap = columns[col].pref_samples_end.at(run_idx - 1);
          auto rank = this->phi->phi_inv_vec[above_old_top_hap].rank1(col);
          this->phi->phi_inv_supp[above_old_top_hap].set(rank, M);

          // set the new sample beg
          columns[col].pref_samples_beg.set(run_idx, M);
          this->phi->phi_vec[M].set(col, true);
          temp_supp.push_back(above_old_top_hap);
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

      // find hapID that is before and after inserted haplotype
      unsigned int hap_before = this->phi->phi(hap_id, col).value();
      unsigned int hap_after = hap_id;

      // Insert and set bit vecs for for hap_before
      auto rank_hap_before = this->phi->phi_inv_vec[hap_before].rank1(col);
      this->phi->phi_inv_supp[hap_before].insert(rank_hap_before, M);
      this->phi->phi_inv_vec[hap_before].set(col, true);

      // Insert and set bit vecs for for hap_after
      auto rank_hap_after = this->phi->phi_vec[hap_after].rank1(col);
      this->phi->phi_supp[hap_after].insert(rank_hap_after, M);
      this->phi->phi_vec[hap_after].set(col, true);

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
      // this->phi->phi_supp[M].push_back(hap_before);
      // this->phi->phi_inv_supp[M].push_back(hap_after);
      temp_supp.push_back(hap_before);
      temp_inv_supp.push_back(hap_after);

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
      vector<pair<unsigned int, unsigned int>> insertion_indices; // stores {index, hapid}
      insertion_indices.emplace_back(M, M);

      // calculate insertion indices
      for(unsigned int col = 1; col < query.size(); ++col) {
        auto retval = w_mod(insertion_indices[col-1].first, col-1,  query[col-1], insertion_indices[col-1].second);
        // cout << retval.first << " : " << retval.second << "\n";
        insertion_indices.push_back(retval);
      }

      // initialize new phi structure
      suc_bv tmp_b, tmp_e;
      for(unsigned int col = 0; col < query.size(); ++col) {
        tmp_b.push_back(false);
        tmp_e.push_back(false);
      }
      this->phi->phi_vec.push_back(tmp_b);
      this->phi->phi_inv_vec.push_back(tmp_e);
      packed_spsi temp_supp;
      packed_spsi temp_inv_supp;

      // perform insertion
      for(unsigned int col = 0; col < query.size(); ++col) {
        Insert(col, insertion_indices[col].first, insertion_indices[col].second, query[col], temp_supp, temp_inv_supp);
      }

      // handling for col == N case
      temp_supp.push_back(this->phi->phi(insertion_indices[N-1].second, N-1).value());
      temp_inv_supp.push_back(insertion_indices[N-1].second);

      this->phi->phi_supp.push_back(temp_supp);
      this->phi->phi_inv_supp.push_back(temp_inv_supp);
      // increment total # of haplotypes
      ++M;
    }

    // Build the reference panel
    void BuildFromVCF(std::string& filename, bool verbose) {
      std::string line= "##";
      if (std::ifstream inFile(filename); inFile.is_open()){
        // go through headers
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
          // M += 2;
          // individual_ids.push_back(token);
        }

        // go through all sites
        vector<int> u, v;
        vector<int> freq;
        unsigned int total_runs = 0;
        int col = 0;
        int cnt = 1;
        bool prev_allele = false;
        // TODO: Possible optimization using sdsl int vec
        vector<int> prefix_arr(M, 0);
        std::iota(prefix_arr.begin(), prefix_arr.end(), 0);
        // TODO: Possible optimization using sdsl int vec
        vector<vector<unsigned int>> sites_where_sample_beg(M);
        vector<vector<unsigned int>> sites_where_sample_end(M);

        while(getline(inFile, line)){
          std::istringstream iss(line);
          token = "";

          // read through the line and store in single_col
          std::vector<bool> single_col;
          for(int i = 0; i < (this->M/2) +9 ; ++i){
            iss >> token;
            if (i < 9){
              continue;
            }
            single_col.push_back(static_cast<bool>(token[0] - '0'));
            single_col.push_back(static_cast<bool>(token[2] - '0'));
          }

          packed_spsi temp_zeros;
          packed_spsi temp_ones;
          packed_spsi temp_combined;
          packed_spsi temp_sample_beg;
          packed_spsi temp_sample_end;
          bool start_with_zero = false;

          assert(single_col.size() == this->M);
          for (int i = 0; i < this->M; ++i) {
            // first allele
            if (i == 0) {
              if (single_col[prefix_arr[i]]) { // allele: 1
                v.push_back(prefix_arr[i]);
              } else { // allele: 0
                u.push_back(prefix_arr[i]);
                start_with_zero = true;
              }
              temp_sample_beg.push_back(prefix_arr[i]);
              prev_allele = single_col[prefix_arr[i]];
              cnt = 1;
              continue;
            }

            // at run-change
            if (single_col[prefix_arr[i]] != prev_allele) {
              freq.push_back(cnt);
              prev_allele = single_col[prefix_arr[i]];
              cnt = 1;
              temp_sample_beg.push_back(prefix_arr[i]);
              temp_sample_end.push_back(prefix_arr[i - 1]);
            } else {
              ++cnt;
            }

            if (single_col[prefix_arr[i]]) { // allele 1
              v.push_back(prefix_arr[i]);
            } else { // allele 0
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
        }
        assert(columns.size() == N);
        // build phi data-structure
        this->phi = new phi_ds(columns, M, N, sites_where_sample_beg, sites_where_sample_end, prefix_arr, verbose);
        assert(col == N);
        total_runs += freq.size();
        cout << "Phi support size (in bytes) = " << this->phi->size_in_bytes(verbose) << "\n";
        cout << "Avg runs = " << static_cast<float>(total_runs)/N << "\n";

        inFile.close();
      }else{
        std::cerr << "Couldn't find : " << filename << "\n";
        exit(1);
      }

    }

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
