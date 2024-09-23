//
// Created by shakyap on 8/20/24.
//

#ifndef DCPBWT_INCLUDE_UTILS_H_
#define DCPBWT_INCLUDE_UTILS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "htslib/vcf.h"

using namespace std;


inline pair<unsigned int, unsigned int> ReadVCF(string &filename) {
  unsigned int M = 0;
  unsigned int N = 0;
  std::string line = "##";
  std::ifstream inFile(filename);
  if (inFile.is_open()) {
    while (line[1] == '#') {
      getline(inFile, line);
    }

    // read individual ids
    std::istringstream iss(line);
    std::string token;
    for (int i = 0; i < 9; ++i) {
      iss >> token;
    }
    while (iss >> token) {
      M += 2;
      // individual_ids.push_back(token);
    }

    // go through all sites
    while (getline(inFile, line)) {
      ++N;
    }
    inFile.close();
  } else {
    std::cerr << "Couldn't find : " << filename << "\n";
    exit(1);
  }
  std::cout << "M (# of haplotypes) = " << M << " : N (# of sites) = " << N << "\n";
  return {M, N};
}

inline void ReadQueryFile(const char* query_vcf, std::vector<std::vector<bool>> &queries){
  htsFile *fp = hts_open(query_vcf, "rb");
  std::cout << "Reading VCF file...\n";
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();
  bool first = true;
  while (bcf_read(fp, hdr, rec) >= 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    // read SAMPLE
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    int max_ploidy = ngt / nsmpl;
    for (i = 0; i < nsmpl; i++) {
      int32_t *ptr = gt_arr + i * max_ploidy;
      for (j = 0; j < max_ploidy; j++) {
        // if true, the sample has smaller ploidy
        if (ptr[j] == bcf_int32_vector_end) break;
        // missing allele
        if (bcf_gt_is_missing(ptr[j])) exit(-1);
        // the VCF 0-based allele index
        int allele_index = bcf_gt_allele(ptr[j]);
        if (first){
          queries.push_back({allele_index == 1});
        } else {
          queries[2*i + j].push_back(allele_index == 1);
        }
      }
    }
    first = false;
    free(gt_arr);
  }
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  bcf_destroy(rec);
}

inline void ReadQueryVCF(string &filename, vector<vector<bool>> &alleles) {
  int M = 0;
  int N = 0;
  std::string line = "##";
  std::ifstream inFile(filename);
  bool first = true;
  if (inFile.is_open()) {
    while (line[1] == '#') {
      getline(inFile, line);
    }

    // read individual ids
    std::istringstream iss(line);
    std::string token;
    for (int i = 0; i < 9; ++i) {
      iss >> token;
    }
    while (iss >> token) {
      M += 2;
      // individual_ids.push_back(token);
    }

    // go through all sites
    while (getline(inFile, line)) {
      ++N;
      std::istringstream iss(line);
      token = "";
//      std::vector<bool> single_col;
      int hap_idx = 0;
      for (int i = 0; i < (M / 2) + 9; ++i) {
        iss >> token;
        if (i < 9) {
          continue;
        }
        if (first){
          alleles.push_back({static_cast<bool>(token[0] - '0')});
          alleles.push_back({static_cast<bool>(token[2] - '0')});
          continue;
        }
        alleles[2*hap_idx].push_back(static_cast<bool>(token[0] - '0'));
        alleles[2*hap_idx+ 1].push_back(static_cast<bool>(token[2] - '0'));
        ++hap_idx;
//        single_col.push_back(static_cast<bool>(token[0] - '0'));
//        single_col.push_back(static_cast<bool>(token[2] - '0'));
      }
      first = false;
//      alleles.push_back(single_col);
    }
    inFile.close();
  } else {
    std::cerr << "Couldn't find : " << filename << "\n";
    exit(1);
  }
//  std::cout << "M (# of haplotypes) = " << M << " : N (# of sites) = " << N << "\n";
}

inline void ReadVCF(string &filename, vector<vector<bool>> &alleles) {
  int M = 0;
  int N = 0;
  std::string line = "##";
  std::ifstream inFile(filename);
  if (inFile.is_open()) {
    while (line[1] == '#') {
      getline(inFile, line);
    }

    // read individual ids
    std::istringstream iss(line);
    std::string token;
    for (int i = 0; i < 9; ++i) {
      iss >> token;
    }
    while (iss >> token) {
      M += 2;
      // individual_ids.push_back(token);
    }

    // go through all sites
    while (getline(inFile, line)) {
      ++N;
      std::istringstream iss(line);
      token = "";
      std::vector<bool> single_col;
      for (int i = 0; i < (M / 2) + 9; ++i) {
        iss >> token;
        if (i < 9) {
          continue;
        }
        single_col.push_back(static_cast<bool>(token[0] - '0'));
        single_col.push_back(static_cast<bool>(token[2] - '0'));
      }
      alleles.push_back(single_col);
    }
    inFile.close();
  } else {
    std::cerr << "Couldn't find : " << filename << "\n";
    exit(1);
  }
  std::cout << "M (# of haplotypes) = " << M << " : N (# of sites) = " << N << "\n";
}



#endif //DCPBWT_INCLUDE_UTILS_H_
