//
// Created by shakyap on 8/9/24.
//

#ifndef PHI_H
#define PHI_H

#include <iostream>
#include <unordered_map>
#include <optional>

#include "dcpbwt_column.h"
#include "dynamic/dynamic.hpp"

using namespace dyn;

class phi_ds {
 public:
  unsigned int total_haplotypes{};
  unsigned int total_sites{};
  /**
   * dynamic sparse bitvectors for phi function
   */
  std::vector<suc_bv> phi_vec;

  /**
   * dynamic sparse bitvectors for phi_inv function
   */
  std::vector<suc_bv> phi_inv_vec;

  /**
   * int vector for prefix samples used by phi function
   */
  std::vector<packed_spsi> phi_supp;

  /**
   * int vector for prefix samples used by phi_inv function
   */
  std::vector<packed_spsi> phi_inv_supp;
  std::vector<packed_spsi> phi_supp_lcp;

  phi_ds() = default;
  ~phi_ds()=default;

  /*
   * constructor of the phi/phi_inv support data structure
   */
  phi_ds(vector<dcpbwt_column> &columns, unsigned int M,
         unsigned int N,
         std::vector<std::vector<unsigned int>> &site_where_samples_beg,
         std::vector<std::vector<unsigned int>> &site_where_samples_end,
         vector<unsigned int> &last_pref,
         vector<unsigned int> &last_div) {
    this->total_haplotypes = M;
    this->total_sites = N;
    this->phi_vec = std::vector<suc_bv>(total_haplotypes);
    this->phi_inv_vec = std::vector<suc_bv>(total_haplotypes);
    this->phi_supp = std::vector<packed_spsi>(total_haplotypes);
    this->phi_inv_supp = std::vector<packed_spsi>(total_haplotypes);
    this->phi_supp_lcp = std::vector<packed_spsi>(total_haplotypes);

    // build sparse bitvector
    for (unsigned int hap = 0; hap < total_haplotypes; hap++) {
      suc_bv tmp_begin, tmp_end;
      // initialize both suc_bv's with 0s (i.e. false)
      for (unsigned int i = 0; i < total_sites + 1; ++i) {
        tmp_begin.push_back(false);
        tmp_end.push_back(false);
      }

      // set the site where the haplotype was at the beginning (or end) of a run
      for (unsigned int k = 0; k < site_where_samples_beg[hap].size(); k++) {
        tmp_begin.set(site_where_samples_beg[hap][k], true);
      }
      for (unsigned int k = 0; k < site_where_samples_end[hap].size(); k++) {
        tmp_end.set(site_where_samples_end[hap][k], true);
      }

      // phi_vec.push_back(tmp_begin);
      phi_vec[hap] = tmp_begin;
      // phi_inv_vec.push_back(tmp_end);
      phi_inv_vec[hap] = tmp_end;
    }

    // iterate over columns
    for (unsigned int col = 0; col < columns.size(); col++) {
      for (unsigned int j = 0; j < columns[col].pref_samples_beg.size(); j++) {
        // Use sample beg to compute phi panel and,
        // use sample_end to compute phi support panel.
        // If are at the top of a column, the haplotype above it doesn't exist.
        // So, use the same haplotype to indicate this case.
        if (j == 0) {
          this->phi_supp[columns[col].pref_samples_beg[j]].push_back(columns[col].pref_samples_beg[j]);
          this->phi_supp_lcp[columns[col].pref_samples_beg[j]].push_back(col);
        } else {
          this->phi_supp[columns[col].pref_samples_beg[j]].push_back(columns[col].pref_samples_end[j - 1]);
          this->phi_supp_lcp[columns[col].pref_samples_beg[j]].push_back(
                   columns[col].div_samples_beg[j]);
        }
        // Use sample end to compute phi_inv panel, and
        // use sample_beg to compute phi_inv support panel
        // If at the bottom of a column, the haplotype below it doesn't exist.
        // So, use the same haplotype to indicate this case.
        if (j == columns[col].pref_samples_beg.size() - 1) {
          this->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
            columns[col].pref_samples_end[j]);
        } else {
          this->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
            columns[col].pref_samples_beg[j + 1]);
        }
      }
    }

    // Use the last prefix array to compute the remaining values for the phi support data structure
    for (unsigned int j = 0; j < total_haplotypes; j++) {
      // Update phi_supp
      if (j == 0) {
        if ((this->phi_supp[last_pref[j]].size() == 0) ||
          this->phi_supp[last_pref[j]].at(this->phi_supp[last_pref[j]].size() - 1)
            != last_pref[j]) { // checks if it's already at the top prior to the last column
          this->phi_supp[last_pref[j]].push_back(last_pref[j]);
          this->phi_supp_lcp[last_pref[j]].push_back(total_sites);
        }
      } else {
        if ((this->phi_supp[last_pref[j]].size() == 0) ||
          this->phi_supp[last_pref[j]].at(this->phi_supp[last_pref[j]].size() - 1) != last_pref[j - 1]) {
          this->phi_supp[last_pref[j]].push_back(last_pref[j - 1]);
          this->phi_supp_lcp[last_pref[j]].push_back(last_div[j]);
        }
      }

      // Update phi_inv_supp
      if (j == this->phi_supp.size() - 1) {
        if (this->phi_inv_supp[last_pref[j]].size() == 0 ||
          this->phi_inv_supp[last_pref[j]].at(this->phi_inv_supp[last_pref[j]].size() - 1) != last_pref[j]) {
          this->phi_inv_supp[last_pref[j]].push_back(last_pref[j]);
        }
      } else {
        if (this->phi_inv_supp[last_pref[j]].size() == 0 ||
          this->phi_inv_supp[last_pref[j]].at(this->phi_inv_supp[last_pref[j]].size() - 1) != last_pref[j + 1]) {
          this->phi_inv_supp[last_pref[j]].push_back(last_pref[j + 1]);
        }
      }
      assert(phi_inv_supp[last_pref[j]].size() > 0);
    }
  }

  std::optional<unsigned int> phi(unsigned int pref, unsigned int col) {
    // assert(col < this->total_sites);
    auto tmp_col = this->phi_vec[pref].rank1(col);
    if (tmp_col == this->phi_supp[pref].size()) {
      tmp_col--;
    }
    auto res = static_cast<unsigned int>(this->phi_supp[pref].at(tmp_col));
    return res;
  }

  std::optional<unsigned int> phi_inv(unsigned int pref, unsigned int col) {
    auto tmp_col = this->phi_inv_vec[pref].rank1(col);
    if (tmp_col == this->phi_inv_supp[pref].size()) {
      tmp_col--;
    }
    auto res = static_cast<unsigned int>(this->phi_inv_supp[pref][tmp_col]);
    return res;
  }

   unsigned int plcp(unsigned int pref, unsigned int col) {
       if (col == 0 || pref == UINT_MAX || (phi(pref, col).value() == phi_inv(pref, col).value()) || (phi(pref, col).value() == pref)) {
           return col;
       }
       auto tmp_col = this->phi_vec[pref].rank1(col);
       if(tmp_col == this->phi_supp[pref].size()){
           tmp_col--;
       }
       unsigned int match_start = 0;
       if (tmp_col > 0){
         match_start = static_cast<int>(this->phi_supp_lcp[pref].at(tmp_col));

       } else{
         assert(tmp_col == 0);
         match_start = static_cast<int>(this->phi_supp_lcp[pref].at(tmp_col));
       }
      return match_start;
   }

  /**
   * function to obtain size in bytes of the phi/phi_inv support data
  */
  unsigned long long size_in_bytes(bool verbose = false) {
    unsigned long long size = 0;
    for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
      size += phi_vec[i].bit_size();
      size += phi_inv_vec[i].bit_size();
      size += phi_supp[i].bit_size();
      size += phi_inv_supp[i].bit_size();
      size += phi_supp_lcp[i].bit_size();
    }
    // convert to bytes
    size /= 8;
    size += sizeof(unsigned int);
    size += sizeof(unsigned int);
    if (verbose) {
      std::cout << "phi support: " << size << " bytes\n";
    }
    return size;
  }
};

#endif //PHI_H
