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
