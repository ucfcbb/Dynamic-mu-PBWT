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
#include <deque>

#include "utils.h"

using namespace dyn;
using namespace std;

//typedef spsi<packed_vector, 256, 2> my_spsi;



class DCPBWT {
 public:
  unsigned int M; // #haplotypes
  unsigned int N; // #sites
  set<unsigned int> haplotype_ids;
  std::vector<dcpbwt_column> columns;
  phi_ds *phi = nullptr;

  // constructor
  DCPBWT(std::string ref_vcf_file, bool verbose) {
    // std::vector<std::vector<bool>> alleles;
    // extract alleles from VCF
    auto retval = ReadVCF(ref_vcf_file);
    this->M = retval.first;
    this->N = retval.second;
    for (unsigned int i = 0; i < M; ++i) {
      haplotype_ids.insert(i);
    }
//    cout << haplotype_ids.size() << "\n";

    // build the ref panel
    // Build(alleles);
    BuildFromVCF(ref_vcf_file, verbose);
  }
  ~DCPBWT() {
    delete phi;
  }

  [[nodiscard]] unsigned int get_run_idx(const unsigned int col, const unsigned int i) const {
    if (i >= (this->M - 1)) {
      return columns[col].combined.size() - 1;
    }
    return columns[col].combined.search(i + 1);
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
//    if (run_idx % 2 == 0) {
    if (!(run_idx & 1)) {
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

  std::pair<unsigned int, unsigned int> column_end_lf(unsigned int col, bool val) {
    if (val) { return {this->M, UINT_MAX}; } // last row, invalid hapID
    unsigned int ind = columns[col].num_zeros;
    unsigned int ridx = 0; // assume column starts with a run of ones
    if (columns[col].start_with_zero) { // the second run at index 1 will be run of ones
      ridx = 1;
    }
    return {ind, columns[col].pref_samples_beg.at(ridx)};
  }

  std::pair<int, int> w_mod(const unsigned int i, const unsigned int col, const bool val, const unsigned int pref_val) {
    if (i == this->M) { return column_end_lf(col, val); }
    const unsigned int run_idx = get_run_idx(col, i);

    if (const bool rval = get_run_val(col, run_idx); rval == val) {
      return {lf(col, i, val), pref_val};
    }

    if (run_idx + 1 >= columns[col].combined.size())
      return column_end_lf(col, val);
    return {lf(col, get_run_head(col, run_idx + 1), val),
            columns[col].pref_samples_beg[run_idx + 1]}; // sample_beg fetch pref-val
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

    auto uv = uv_trick(col, i);

    if (value) {
      return columns[col].num_zeros + uv.second + offset;
    }
    return uv.first + offset;
  }

  unsigned int reverse_lf(const unsigned int col, unsigned int idx) {
    if (col == 0) {
      return 0;
    }
    unsigned int prev_col = col - 1;
    unsigned int num_zeros = this->columns[prev_col].num_zeros;
    // came from a run of zeros in prev col
//    bool found = false;
    unsigned offset = 0;
    if (idx < num_zeros) {
      unsigned int zero_idx = this->columns[prev_col].zeros.search(idx + 1);
      unsigned int zeros_before = 0;
      if (zero_idx > 0)
        zeros_before = this->columns[prev_col].zeros.psum(zero_idx - 1);
      offset = idx - zeros_before;
      unsigned int run_idx = 0;
      if (this->columns[prev_col].start_with_zero) {
        run_idx = zero_idx * 2;
      } else {
        run_idx = zero_idx * 2 + 1;
      }
      unsigned int run_head = 0;
      if (run_idx > 0)
        run_head = this->columns[prev_col].combined.psum(run_idx - 1);
      return run_head + offset;
    } else {
      idx -= num_zeros;

      unsigned int ones_idx = this->columns[prev_col].ones.search(idx + 1);
      unsigned int run_idx = 0;
      if (this->columns[prev_col].start_with_zero) {
        run_idx = ones_idx * 2 + 1;
      } else {
        run_idx = ones_idx * 2;
      }

      unsigned int ones_before = 0;
      if (ones_idx > 0)
        ones_before = this->columns[prev_col].ones.psum(ones_idx - 1);
      offset = idx - ones_before;
      unsigned int run_head = 0;
      if (run_idx > 0)
        run_head = this->columns[prev_col].combined.psum(run_idx - 1);
      return run_head + offset;
    }
  }

  unsigned int zeros_before(const unsigned int col, const unsigned int run_idx) {
    unsigned int retval = 0;
    assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
    if (columns[col].start_with_zero) {
      if (run_idx % 2 == 0)
        retval = columns[col].zeros.psum((run_idx - 2) / 2);
      else
        retval = columns[col].zeros.psum((run_idx - 1) / 2);
    } else {
      if (run_idx % 2 == 0)
        retval = columns[col].zeros.psum((run_idx - 1) / 2);
      else
        retval = columns[col].zeros.psum((run_idx - 2) / 2);
    }
    return retval;
  }

  unsigned int ones_before(const unsigned int col, const unsigned int run_idx) {
    unsigned int retval = 0;
    assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
    if (columns[col].start_with_zero) {
      if (run_idx % 2 == 0)
        retval = columns[col].ones.psum((run_idx - 1) / 2);
      else
        retval = columns[col].ones.psum((run_idx - 2) / 2);
    } else {
      if (run_idx % 2 == 0)
        retval = columns[col].ones.psum((run_idx - 2) / 2);
      else
        retval = columns[col].ones.psum((run_idx - 1) / 2);
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
      if (columns[col].start_with_zero) {
        return make_pair(zeros_before(col, run_idx), 0);
      }
      return make_pair(0, ones_before(col, run_idx));
    }
    u = zeros_before(col, run_idx);
    v = ones_before(col, run_idx);
    return make_pair(u, v);
  }

  [[nodiscard]] bool AlleleMatchesRun(const unsigned int col, const unsigned int run_idx, const bool allele) const {
    bool match = false;
    if (this->columns[col].start_with_zero) {
      // even run & allele has to be zero
      if (run_idx & 1) { // odd run
        if (allele)
          match = true;
      } else {
        if (!allele)
          match = true;
      }
    } else {
      // even run & allele has to be one
      if (run_idx & 1) { // odd run
        if (!allele) {
          match = true;
        }
      } else { // even run
        if (allele)
          match = true;
      }
    }
    return match;
  }

  void InsertAtRunStart(const unsigned int col,
                        const unsigned idx,
                        const unsigned int hap_id,
                        bool allele,
                        packed_spsi &temp_supp,
                        packed_spsi &temp_inv_supp,
                        packed_spsi &temp_div_supp,
                        vector<unsigned int> &temp_div_query,
                        vector<unsigned int> &temp_div_below_query) {

    unsigned int run_idx = get_run_idx(col, idx);
    unsigned int hap_before = 0;
    if (idx == 0) {
      hap_before = UINT_MAX;
    } else {
      hap_before = this->columns[col].pref_samples_end[run_idx - 1];
    }

    if (AlleleMatchesRun(col, run_idx, allele)) {
      // update run info
      this->columns[col].combined.increment(run_idx, 1);
      if (allele) {
        this->columns[col].ones.increment(run_idx / 2, 1);
      } else {
        this->columns[col].zeros.increment(run_idx / 2, 1);
        ++this->columns[col].num_zeros;
      }

      // update new run head
      this->columns[col].pref_samples_beg.set(run_idx, this->M);

      // unset old head of run
      auto col_rank = this->phi->phi_vec[hap_id].rank1(col);
      this->phi->phi_supp[hap_id].remove(col_rank);
      this->phi->phi_vec[hap_id].set(col, false);

      // set phi for new run head
      this->phi->phi_vec[this->M].set(col, true);
      temp_supp.push_back(hap_before);

      // update phi for hap_before as well
      if (hap_before != UINT_MAX) {
        col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
        this->phi->phi_inv_supp[hap_before].set(col_rank, this->M);
      }

      // TODO: Update div sample
      // Replaces the previous div sample
      this->columns[col].div_samples_beg.set(run_idx, temp_div_query[col]);
      temp_div_supp.push_back(temp_div_query[col]);
    } else {
      // If at the very top of a column
      if (idx == 0) {
        // Create a new Run
        // Insert new run info
        assert(run_idx == 0);
        this->columns[col].combined.insert(run_idx, 1);
        if (allele) {
          this->columns[col].ones.insert(run_idx, 1);
          this->columns[col].start_with_zero = false;
        } else {
          this->columns[col].zeros.insert(run_idx, 1);
          this->columns[col].start_with_zero = true;
          ++this->columns[col].num_zeros;
        }

        // insert into pref_beg
        this->columns[col].pref_samples_beg.insert(run_idx, this->M);
        this->columns[col].pref_samples_end.insert(run_idx, this->M);

        // update phi for old run head
        auto col_rank = this->phi->phi_vec[hap_id].rank1(col);
        this->phi->phi_supp[hap_id].set(col_rank, this->M);

        // update for newly inserted head
        this->phi->phi_vec[this->M].set(col, true);
        this->phi->phi_inv_vec[this->M].set(col, true);
        temp_supp.push_back(hap_before);
        temp_inv_supp.push_back(hap_id);

        // TODO: Update div sample
        // Insert new div sample and also update the one that's below it
        // Update the one that's below first
        this->columns[col].div_samples_beg.set(run_idx, temp_div_below_query[col]);
        this->columns[col].div_samples_beg.insert(0, temp_div_query[col]);
        temp_div_supp.push_back(col);
      } else {
        /* No need to create a new run
         * Simply insert into previous run
         */
        // update run info
        this->columns[col].combined.increment(run_idx - 1, 1);
        if (allele) {
          this->columns[col].ones.increment((run_idx - 1) / 2, 1);
        } else {
          this->columns[col].zeros.increment((run_idx - 1) / 2, 1);
          ++this->columns[col].num_zeros;
        }

        // Update pref_end (i.e. inserted hap will be at the end of prev run)
        this->columns[col].pref_samples_end.set(run_idx - 1, this->M);

        // update phi for old hap at prev run's end
        auto col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
        this->phi->phi_inv_supp[hap_before].remove(col_rank);
        this->phi->phi_inv_vec[hap_before].set(col, false);

        // update for newly inserted hap
        this->phi->phi_inv_vec[this->M].set(col, true);
        temp_inv_supp.push_back(hap_id);

        // update phi_inv value for hap_id
        if (hap_id != UINT_MAX) {
          col_rank = this->phi->phi_vec[hap_id].rank1(col);
          this->phi->phi_supp[hap_id].set(col_rank, this->M);

          // TODO: Update div sample of the seq below
          this->columns[col].div_samples_beg.set(run_idx, temp_div_below_query[col]);
        }
      }
    }
  }

  void InsertAtBottom(const unsigned int col,
                      const unsigned int run_idx,
                      const unsigned int hap_id,
                      const bool allele,
                      packed_spsi &temp_supp,
                      packed_spsi &temp_inv_supp,
                      packed_spsi &temp_div_supp,
                      vector<unsigned int> &temp_div_query,
                      vector<unsigned int> &temp_div_below_query) {
    if (AlleleMatchesRun(col, run_idx, allele)) {
      // Update Run info
      columns[col].combined.increment(run_idx, 1);
      if (allele) {
        columns[col].ones.increment(run_idx / 2, 1);
      } else {
        columns[col].zeros.increment(run_idx / 2, 1);
        ++columns[col].num_zeros;
      }

      // Update pref end
      unsigned int old_bottom = this->columns[col].pref_samples_end.at(run_idx);
      this->columns[col].pref_samples_end.set(run_idx, this->M);

      // Update old bottom's Phi structures
      auto col_rank = this->phi->phi_inv_vec[old_bottom].rank1(col);
      this->phi->phi_inv_supp[old_bottom].remove(col_rank);
      this->phi->phi_inv_vec[old_bottom].set(col, false);

      // Update new bottom's Phi structures
      this->phi->phi_inv_vec[this->M].set(col, true);
      temp_inv_supp.push_back(UINT_MAX);
    } else {
      // Insert run info
      this->columns[col].combined.push_back(1);
      if (allele) {
        this->columns[col].ones.push_back(1);
      } else {
        this->columns[col].zeros.push_back(1);
        ++this->columns[col].num_zeros;
      }

      assert(run_idx == this->columns[col].pref_samples_end.size() - 1);
      unsigned int hap_before = this->columns[col].pref_samples_end.at(run_idx);

      // update pref beg/end
      this->columns[col].pref_samples_beg.push_back(this->M);
      this->columns[col].pref_samples_end.push_back(this->M);

      // Update phi supports for inserted haplotype
      this->phi->phi_vec[this->M].set(col, true);
      this->phi->phi_inv_vec[this->M].set(col, true);
      temp_supp.push_back(hap_before);
      temp_inv_supp.push_back(UINT_MAX);

      // Update for hap_before
      auto col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
      this->phi->phi_inv_supp[hap_before].set(col_rank, this->M);

      // TODO: Insert Div val since this is a new run
      this->columns[col].div_samples_beg.push_back(temp_div_query[col]);
      temp_div_supp.push_back(temp_div_query[col]);
    }
  }

  void Insert(const unsigned int col,
              const unsigned idx,
              const unsigned int hap_id,
              bool allele,
              packed_spsi &temp_supp,
              packed_spsi &temp_inv_supp,
              packed_spsi &temp_div_supp,
              vector<unsigned int> &temp_div_query,
              vector<unsigned int> &temp_div_below_query) {
    // assumes Non-empty ref panel (i.e. a panel is built prior to insertion)
    auto run_idx = get_run_idx(col, idx);
    if (isRunStart(col, run_idx, hap_id)) {
      InsertAtRunStart(col, idx, hap_id, allele, temp_supp, temp_inv_supp, temp_div_supp, temp_div_query, temp_div_below_query);
      return;
    }

    if (idx == this->M) {
      InsertAtBottom(col, run_idx, hap_id, allele, temp_supp, temp_inv_supp, temp_div_supp, temp_div_query, temp_div_below_query);
      return;
    }

    // if run value is same as the allele
    if (AlleleMatchesRun(col, run_idx, allele)) {
      columns[col].combined.increment(run_idx, 1);
      if (allele) {
        columns[col].ones.increment(run_idx / 2, 1);
      } else {
        columns[col].zeros.increment(run_idx / 2, 1);
        ++columns[col].num_zeros;
      }
    } else {
      /*
       * Create a new run as the existing run will be split
       */
      auto hap_before_opt = this->phi->phi(hap_id, col);
      assert(hap_before_opt.has_value());
      unsigned int hap_before = hap_before_opt.value();
      unsigned int hap_after = hap_id;

      // Update run info
      unsigned int run_head = get_run_head(col, run_idx);
      unsigned int left_rem = idx - run_head;
      unsigned int right_rem = columns[col].combined.at(run_idx) - left_rem;
      columns[col].combined.set(run_idx, left_rem);
      columns[col].combined.insert(run_idx + 1, 1);
      columns[col].combined.insert(run_idx + 2, right_rem);

      if (allele) {
        // Inserting 1 in a run of 0s
        if (run_idx == 0) {
          this->columns[col].ones.insert(run_idx, 1);
        } else {
//          assert(run_idx > 0);
          this->columns[col].ones.insert((run_idx - 1) / 2 + 1, 1);
        }

        // Update Zeros vector
        this->columns[col].zeros.set(run_idx / 2, left_rem);
        if (run_idx / 2 + 1 >= this->columns[col].zeros.size())
          this->columns[col].zeros.push_back(right_rem);
        else
          this->columns[col].zeros.insert(run_idx / 2 + 1, right_rem);
      } else {
        // Inserting 0 in a run of 1s
        if (run_idx == 0) {
          this->columns[col].zeros.insert(run_idx, 1);
        } else {
          this->columns[col].zeros.insert((run_idx - 1) / 2 + 1, 1);
        }
        ++this->columns[col].num_zeros;

        // Update Ones vector
        this->columns[col].ones.set(run_idx / 2, left_rem);
        if (run_idx / 2 + 1 >= this->columns[col].ones.size())
          this->columns[col].ones.push_back(right_rem);
        else
          this->columns[col].ones.insert(run_idx / 2 + 1, right_rem);
      }

      // Update pref_beg/end
      // Update Pref Beg
      if (run_idx + 1 >= this->columns[col].pref_samples_beg.size()) {
        this->columns[col].pref_samples_beg.push_back(this->M);
        this->columns[col].pref_samples_beg.push_back(hap_after);
      } else {
        this->columns[col].pref_samples_beg.insert(run_idx + 1, this->M);
        this->columns[col].pref_samples_beg.insert(run_idx + 2, hap_after);
      }
      // Update Pref end
      this->columns[col].pref_samples_end.insert(run_idx, this->M);
      this->columns[col].pref_samples_end.insert(run_idx, hap_before);

      // Update phi supports for hap_before
      auto col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
      this->phi->phi_inv_supp[hap_before].insert(col_rank, this->M);
      this->phi->phi_inv_vec[hap_before].set(col, true);

      // Update phi supports for hap_after
      col_rank = this->phi->phi_vec[hap_after].rank1(col);
      this->phi->phi_supp[hap_after].insert(col_rank, this->M);
      this->phi->phi_vec[hap_after].set(col, true);
      this->phi->phi_supp[hap_after].insert(col_rank, temp_div_below_query[col]);

      // update phi structure for Inserted Haplotype
      this->phi->phi_vec[M].set(col, true);
      this->phi->phi_inv_vec[M].set(col, true);
      temp_supp.push_back(hap_before);
      temp_inv_supp.push_back(hap_after);

      // TODO: Insert div value for inserted hap and the hap below it since we're forming two new head of runs
      if (run_idx + 1 >= this->columns[col].combined.size()) {
        this->columns[col].div_samples_beg.push_back(temp_div_query[col]);
        this->columns[col].div_samples_beg.push_back(temp_div_below_query[col]);

        this->phi->phi_supp[hap_after].push_back(temp_div_below_query[col]);
        temp_div_supp.push_back(temp_div_query[col]);
      } else {
        this->columns[col].div_samples_beg.insert(run_idx + 1, temp_div_below_query[col]);
        this->columns[col].div_samples_beg.insert(run_idx + 1, temp_div_query[col]);
        temp_div_supp.push_back(temp_div_query[col]);
      }
    }
  }

  /*
     * Given a position in the PBWT:
     * hapInd in [0,height)
     * col in [1, width]
     * output the value haplotype getPrefix[col] has at position col - 1
     */
  bool get_curr_char(unsigned int hapInd, unsigned int col) {
    assert(hapInd >= 0 && hapInd < this->M);
    assert(col != 0 && col <= this->N);
    return hapInd >= this->columns[col - 1].num_zeros;
  }

  void InsertSinglelHaplotype(std::vector<bool> &query) {
    assert(query.size() == N);
    vector<pair<unsigned int, unsigned int>> insertion_indices(this->N + 1); // stores {index, hapid}
    insertion_indices[0].first = this->M;
    insertion_indices[0].second = UINT_MAX;

    // calculate insertion indices
    for (unsigned int col = 0; col < query.size(); ++col) {
      insertion_indices[col + 1] = w_mod(insertion_indices[col].first, col, query[col], insertion_indices[col].second);
    }

    cout << insertion_indices.size() << "\n";
    // TODO: calculate divergence values
    vector<unsigned int> temp_div_query(this->N + 1, 0);
    vector<unsigned int> temp_div_below_query(this->N + 1, 0);
    // update divergence values
    clock_t START2 = clock();
    unsigned int zs = 0; // starting position of match between query and seqn above it
    unsigned int bs = 0; // starting position of match between query and seqn below it
    unsigned int hap_id_above_prev = -1;
    for (unsigned int k = this->N; k > 0; --k) {
      zs = k;
      bs = k;
      unsigned int hap_id = insertion_indices[k].second;
      unsigned int position = insertion_indices[k].first;

      // update divergence values for sequence BELOW query
      // "starting pos" definition of divergence value
      if (k != this->N && insertion_indices[k + 1].second == hap_id)
        bs = temp_div_below_query[k + 1];
      else if (position != this->M) {
        while (bs > 0 &&
          get_curr_char(position, bs) == query[bs - 1]) {
//          position = reverse_lf(bs, position, verbose);
          position = reverse_lf(bs, position);
          --bs;
        }
      }

      // update divergence values for query sequence
      // "starting pos" definition of divergence value
      position = insertion_indices[k].first;
      if (position-- != 0) {
        std::optional<unsigned int> hap_id_above;// = this->phi->phi(hap_id, k);;
        if (position == this->M - 1) {
          hap_id_above = this->columns[k].pref_samples_end.at(this->columns[k].combined.size() - 1);
        } else {
          hap_id_above = this->phi->phi(hap_id, k);
        }
        assert(hap_id_above.has_value());
        if (k != this->N && hap_id_above_prev == hap_id_above.value())
          zs = temp_div_query[k + 1];
        else {
          while (zs > 0 &&
            get_curr_char(position, zs) == query[zs - 1]) {
//            position = reverse_lf(zs, position, verbose);
            position = reverse_lf(zs, position);
            --zs;
          }
          hap_id_above_prev = hap_id_above.value();
        }
      }
      temp_div_query[k] = zs;
      temp_div_below_query[k] = bs;
    }
    temp_div_query[0] = temp_div_below_query[0] = 0;

    // initialize new phi structure
    suc_bv tmp_b, tmp_e;
    for (unsigned int col = 0; col < query.size(); ++col) {
      tmp_b.push_back(false);
      tmp_e.push_back(false);
    }
    this->phi->phi_vec.push_back(tmp_b);
    this->phi->phi_inv_vec.push_back(tmp_e);
    packed_spsi temp_supp;
    packed_spsi temp_inv_supp;
    packed_spsi temp_div_supp;

    // perform insertion
    // where div samples are also updated accordingly
    for (unsigned int col = 0; col < query.size(); ++col) {
      Insert(col,
             insertion_indices[col].first,
             insertion_indices[col].second,
             query[col],
             temp_supp,
             temp_inv_supp,
             temp_div_supp,
             temp_div_query,
             temp_div_below_query);
    }

    // handling for col == N case
    auto hap_above_opt = this->phi->phi(insertion_indices[N - 1].second, N - 1);
    assert(hap_above_opt.has_value());
    auto hap_above = hap_above_opt.value();
    auto hap_below = insertion_indices[N - 1].second;
    temp_supp.push_back(hap_above);
    temp_inv_supp.push_back(hap_below);
    this->phi->phi_supp.push_back(temp_supp);
    this->phi->phi_inv_supp.push_back(temp_inv_supp);
    this->phi->phi_supp_lcp.push_back(temp_div_supp);

    // Update for hap_below and hap_above too
    if (hap_above != UINT_MAX) {
      this->phi->phi_inv_supp[hap_above].set(this->phi->phi_inv_supp[hap_above].size() - 1, this->M);
    }
    if (hap_below != UINT_MAX) {
      this->phi->phi_supp[hap_below].set(this->phi->phi_supp[hap_below].size() - 1, this->M);
    }

    // increment total # of haplotypes
    this->haplotype_ids.insert(this->M);
    ++this->M;
    ++this->phi->total_haplotypes;
  }

  bool isRunWithSingleElement(const unsigned int col, const unsigned int run_idx) {
    return columns[col].combined.at(run_idx) == 1;
  }

  // helper methods
  [[nodiscard]] bool isRunEnd(const unsigned int col, const unsigned int run_idx, const unsigned int hap_id) const {
    return hap_id == this->columns[col].pref_samples_end.at(run_idx);
  }

  [[nodiscard]] bool isRunStart(const unsigned int col, const unsigned int run_idx, const unsigned int hap_id) const {
    return hap_id == this->columns[col].pref_samples_beg.at(run_idx);
  }

  void DeleteSingleRun(const unsigned int col, const unsigned int idx, const unsigned int hap_id, const bool allele) {
    unsigned int run_idx = 0;
//    unsigned int run_idx = get_run_idx(col, idx);
    // if at the top of a column
    if (idx == 0) {
      // new top will be the head of next run
      unsigned int new_top = UINT_MAX;
      if (idx + 1 < this->columns[col].pref_samples_beg.size())
        new_top = this->columns[col].pref_samples_beg.at(idx + 1);

      // update pref_beg/end
      this->columns[col].pref_samples_beg.remove(run_idx);
      this->columns[col].pref_samples_end.remove(run_idx);

      // update run info spsis
      this->columns[col].combined.remove(idx);
      if (allele) {
        this->columns[col].ones.remove(idx);
        this->columns[col].start_with_zero = true;
      } else {
        this->columns[col].zeros.remove(idx);
        --this->columns[col].num_zeros;
        this->columns[col].start_with_zero = false;
      }

      // update phi structures for removed haplotype
      auto col_rank = this->phi->phi_vec[hap_id].rank1(col);
      if (this->phi->phi_supp[hap_id].size() > 1)
        this->phi->phi_supp[hap_id].remove(col_rank);
      col_rank = this->phi->phi_inv_vec[hap_id].rank1(col);
      if (this->phi->phi_inv_supp[hap_id].size() > 1)
        this->phi->phi_inv_supp[hap_id].remove(col_rank);
      this->phi->phi_vec[hap_id].set(col, false);
      this->phi->phi_inv_vec[hap_id].set(col, false);

      // update phi for the new top
      if (new_top != UINT_MAX) {
        col_rank = this->phi->phi_vec[new_top].rank1(col);
        this->phi->phi_supp[new_top].set(col_rank, UINT_MAX);
      }

      // TODO: Update div of the next top haplotype
      // TODO: Update phi_supp_lcp as well
      // remove the current top
      this->columns[col].div_samples_beg.remove(0);
      // Update div of the new top
      this->columns[col].div_samples_beg.set(0, col);

    } else if (idx == this->M - 1) { // at the bottom of a column
      run_idx = this->columns[col].pref_samples_beg.size() - 1;

      // new bottom will be the head of next run
      unsigned int new_bottom = this->columns[col].pref_samples_end[run_idx - 1];

      // update pref_beg/end
      this->columns[col].pref_samples_beg.remove(run_idx);
      this->columns[col].pref_samples_end.remove(run_idx);

      // update run info spsis
      this->columns[col].combined.remove(run_idx);
      if (allele) {
        assert(this->columns[col].ones.at(run_idx / 2) == 1);
        this->columns[col].ones.remove(run_idx / 2);
      } else {
        assert(this->columns[col].zeros.at(run_idx / 2) == 1);
        this->columns[col].zeros.remove(run_idx / 2);
        --this->columns[col].num_zeros;
      }

      // update phi structures for removed haplotype
      auto col_rank = this->phi->phi_vec[hap_id].rank1(col);
      if (this->phi->phi_supp[hap_id].size() > 1)
        this->phi->phi_supp[hap_id].remove(col_rank);
      col_rank = this->phi->phi_inv_vec[hap_id].rank1(col);
      if (this->phi->phi_inv_supp[hap_id].size() > 1)
        this->phi->phi_inv_supp[hap_id].remove(col_rank);
      this->phi->phi_vec[hap_id].set(col, false);
      this->phi->phi_inv_vec[hap_id].set(col, false);

      // update phi for the new bottom
      col_rank = this->phi->phi_inv_vec[new_bottom].rank1(col);
      this->phi->phi_inv_supp[new_bottom].set(col_rank, UINT_MAX);

      // TODO: Remove the last div sample
      // TODO: Update phi_supp_lcp as well
      this->columns[col].div_samples_beg.remove(this->columns[col].div_samples_beg.size() - 1);

    } else {
      run_idx = get_run_idx(col, idx);
      unsigned int hap_before = this->columns[col].pref_samples_end[run_idx - 1];
      unsigned int hap_after = this->columns[col].pref_samples_beg[run_idx + 1];
      unsigned int hap_id_div_val = this->columns[col].div_samples_beg.at(run_idx);
      unsigned int hap_after_div_val = this->columns[col].div_samples_beg.at(run_idx + 1);

      // 0 0 0 0 0 0 0  |1|   0 0 0 0 0 0
      //  ridx-1
      // merge run_idx + 1 into run_idx - 1
      // update run info spsis
      this->columns[col].combined.increment(run_idx - 1, this->columns[col].combined.at(run_idx + 1));
      this->columns[col].combined.remove(run_idx + 1);
      this->columns[col].combined.remove(run_idx);
      if (allele) {
        this->columns[col].zeros.increment((run_idx - 1) / 2, this->columns[col].zeros.at((run_idx + 1) / 2));
        this->columns[col].zeros.remove((run_idx + 1) / 2);

        this->columns[col].ones.remove(run_idx / 2);
      } else {
        this->columns[col].ones.increment((run_idx - 1) / 2, this->columns[col].ones.at((run_idx + 1) / 2));
        this->columns[col].ones.remove((run_idx + 1) / 2);

        this->columns[col].zeros.remove(run_idx / 2);
        --this->columns[col].num_zeros;
      }

      // update pref_beg/end
      this->columns[col].pref_samples_beg.remove(run_idx + 1);
      this->columns[col].pref_samples_beg.remove(run_idx);

      this->columns[col].pref_samples_end.remove(run_idx);
      this->columns[col].pref_samples_end.remove(run_idx - 1);

      // unset phi structures for hap_before and hap_after
      // For hap_before
      auto col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
      this->phi->phi_inv_supp[hap_before].remove(col_rank);
      this->phi->phi_inv_vec[hap_before].set(col, false);
      // For hap_after
      col_rank = this->phi->phi_vec[hap_after].rank1(col);
      this->phi->phi_supp[hap_after].remove(col_rank);
      this->phi->phi_supp_lcp[hap_after].remove(col_rank);
      this->phi->phi_vec[hap_after].set(col, false);

      // unset for hap_id
      col_rank = this->phi->phi_vec[hap_id].rank1(col);
      this->phi->phi_supp[hap_id].remove(col_rank);
      this->phi->phi_supp_lcp[hap_id].remove(col_rank);
      this->phi->phi_vec[hap_id].set(col, false);

      col_rank = this->phi->phi_inv_vec[hap_id].rank1(col);
      this->phi->phi_inv_supp[hap_id].remove(col_rank);
      this->phi->phi_inv_vec[hap_id].set(col, false);

      // TODO: Remove the current run and next run's div samples (Think about how deletion changes the div values)
      // TODO: Update phi_supp_lcp as well
      this->columns[col].div_samples_beg.remove(run_idx);
      this->columns[col].div_samples_beg.remove(run_idx);
    }
  }

  void DeleteAtRunStart(const unsigned int col,
                        const unsigned int run_idx,
                        const unsigned int hap_id,
                        const bool allele) {
    unsigned int hap_before = 0;
    if (run_idx == 0) {
      hap_before = UINT_MAX;
    } else {
      hap_before = this->columns[col].pref_samples_end.at(run_idx - 1);
    }
    auto hap_after_opt = this->phi->phi_inv(hap_id, col);
    assert(hap_after_opt.has_value());
    unsigned int hap_after = hap_after_opt.value();
    auto hap_after_div_val = static_cast<int>(this->phi->plcp(hap_after, col));

    // update run info spsi
    this->columns[col].combined.decrement(run_idx, 1);
    if (allele) {
      this->columns[col].ones.decrement(run_idx / 2, 1);
    } else {
      this->columns[col].zeros.decrement(run_idx / 2, 1);
      --this->columns[col].num_zeros;
    }

    // update pref_beg with new head of run
    this->columns[col].pref_samples_beg.set(run_idx, hap_after);

    // update phi for hap_id
    auto col_rank = this->phi->phi_vec[hap_id].rank1(col);
    this->phi->phi_supp[hap_id].remove(col_rank);
    this->phi->phi_vec[hap_id].set(col, false);

    // update phi for new head of run
    col_rank = this->phi->phi_vec[hap_after].rank1(col);
    this->phi->phi_supp[hap_after].insert(col_rank, hap_before);
    this->phi->phi_vec[hap_after].set(col, true);

    // update phi for hap_before
    if (hap_before != UINT_MAX) {
      col_rank = this->phi->phi_inv_vec[hap_before].rank1(col);
      this->phi->phi_inv_supp[hap_before].set(col_rank, hap_after);
    }

    // TODO: Update div for hap below
    auto curr_div = static_cast<int>(this->columns[col].div_samples_beg.at(run_idx));
    this->columns[col].div_samples_beg.set(run_idx, max(curr_div, hap_after_div_val));
  }

  void DeleteAtRunEnd(const unsigned int col,
                      const unsigned int run_idx,
                      const unsigned int hap_id,
                      const bool allele) {
    auto new_bottom_opt = this->phi->phi(hap_id, col);
    assert(new_bottom_opt.has_value());
    auto new_bottom = new_bottom_opt.value();
    unsigned int below_hap_id = UINT_MAX;
    if (run_idx + 1 < this->columns[col].pref_samples_beg.size()) {
      below_hap_id = this->columns[col].pref_samples_beg.at(run_idx + 1);
    }


    // update pref_end
    this->columns[col].pref_samples_end.set(run_idx, new_bottom);

    // update hap_id
    auto col_rank = this->phi->phi_inv_vec[hap_id].rank1(col);
    this->phi->phi_inv_supp[hap_id].remove(col_rank);
    this->phi->phi_inv_vec[hap_id].set(col, false);

    // update for run_info
    this->columns[col].combined.decrement(run_idx, 1);
    if (allele) {
      this->columns[col].ones.decrement(run_idx / 2, 1);
    } else {
      this->columns[col].zeros.decrement(run_idx / 2, 1);
      --this->columns[col].num_zeros;
    }

    // update for new bottom and the haplotype below hap_id
    col_rank = this->phi->phi_inv_vec[new_bottom].rank1(col);
    this->phi->phi_inv_supp[new_bottom].insert(col_rank, below_hap_id);
    this->phi->phi_inv_vec[new_bottom].set(col, true);

    if (below_hap_id != UINT_MAX) {
      col_rank = this->phi->phi_vec[below_hap_id].rank1(col);
      this->phi->phi_supp[below_hap_id].set(col_rank, new_bottom);

      // TODO: Update div value for the hap below
      auto curr_div_val = this->phi->plcp(hap_id, col);
      if (curr_div_val > this->columns[col].div_samples_beg.at(run_idx + 1)) {
        this->columns[col].div_samples_beg.set(run_idx + 1, curr_div_val);
      }
    }

  }

  void Delete(const unsigned int col, const unsigned int idx, const unsigned int hap_id, const bool allele) {
    unsigned int run_idx = get_run_idx(col, idx);
    // if deleting entire run
    if (isRunWithSingleElement(col, run_idx)) {
      DeleteSingleRun(col, idx, hap_id, allele);
    } else if (isRunStart(col, run_idx, hap_id)) {// if deleting at start of a run
      DeleteAtRunStart(col, run_idx, hap_id, allele);
    } else if (isRunEnd(col, run_idx, hap_id)) {// if deleting at end of a run
      DeleteAtRunEnd(col, run_idx, hap_id, allele);
    } else {
      this->columns[col].combined.decrement(run_idx, 1);
      if (allele) {
        this->columns[col].ones.decrement(run_idx / 2, 1);
      } else {
        this->columns[col].zeros.decrement(run_idx / 2, 1);
        --this->columns[col].num_zeros;
      }

    }
  }

  /*
   * Input: haplotype index which will also be haplotype id for the first column (index is 0-based)
   * Result: Deletes that haplotype's allele from all columns
   */
  void DeleteSingleHaplotype(const unsigned int idx) {
    if (this->M == 0) {
      cerr << "Panel is empty. Nothing to delete! \n";
      return;
    }

    // TODO: need a way to convert the hap_idx to it's corresponding hap_ID if the panel's been updated already (i.e. in case of deletion)
    // E.g. hap_idx: 0 1 2 3 4 5 => delete(3) => 0 1 2 x 3 4
    //       hap_ID: 0 1 2 3 4 5 => delete(3) => 0 1 2 3 4 5
    // Get appropriate haplotype_index
    auto it = haplotype_ids.begin();
    std::advance(it, idx);
    unsigned int hap_id = *it;
//    unsigned int idx = 0;
//    auto ret = find(haplotype_ids.begin(), haplotype_ids.end(), hap_id);
//    idx = distance(haplotype_ids.begin(), ret);

    assert (idx < this->M && idx >= 0);
    vector<pair<unsigned int, unsigned int>> haplotype_info; // {index, hapid, allele}
    vector<bool> alleles;

    // find where each haplotype is mapped in each column
    for (unsigned int col = 0; col < this->N; ++col) {
      if (col == 0) {
        haplotype_info.emplace_back(idx, hap_id);
        continue;
      }
      unsigned int run_idx = get_run_idx(col - 1, haplotype_info[col - 1].first);
      const bool allele = this->get_run_val(col - 1, run_idx);
      auto retval = w_mod(haplotype_info[col - 1].first, col - 1, allele, haplotype_info[col - 1].second);
      haplotype_info.emplace_back(retval);
      alleles.push_back(allele);
    }

    // get last allele
    auto run_idx = get_run_idx(this->N - 1, haplotype_info[this->N - 1].first);
    bool allele = this->get_run_val(this->N - 1, run_idx);
    alleles.push_back(allele);
    assert(alleles.size() == haplotype_info.size());
    assert(alleles.size() == this->N);
    if (this->M == 1) {
      cout << "Only single haplotype remaining to delete!\n";
      for (unsigned int col = 0; col < this->N; ++col) {
        this->columns[col].pref_samples_beg.set(0, 0);
        this->columns[col].pref_samples_end.set(0, 0);
      }
    }
    // delete from each column
    for (unsigned int col = 0; col < this->N; ++col) {
      Delete(col, haplotype_info[col].first, haplotype_info[col].second, alleles[col]);
    }

    // checks
    for (auto i = 0; i < this->N; ++i) {
      assert(static_cast<bool>(this->phi->phi_vec[haplotype_info[0].second].at(i)) == false);
      assert(static_cast<bool>(this->phi->phi_inv_vec[haplotype_info[0].second].at(i)) == false);
    }
    assert(this->phi->phi_supp[haplotype_info[0].second].size() == 1);
    assert(this->phi->phi_inv_supp[haplotype_info[0].second].size() == 1);

    // Fix for last column update
    auto last_above = this->phi->phi_supp[haplotype_info[0].second].at(0);
    auto last_below = this->phi->phi_inv_supp[haplotype_info[0].second].at(0);
    if (last_above != UINT_MAX) {
      this->phi->phi_inv_supp[last_above].set(this->phi->phi_inv_supp[last_above].size() - 1, last_below);
    }
    if (last_below != UINT_MAX) {
      this->phi->phi_supp[last_below].set(this->phi->phi_supp[last_below].size() - 1, last_above);
    }

    // Update phi supports
//    if (idx != this->M - 1) {
//      for (int i = 0; i < this->M; ++i) {
//        if (i == hap_id)
//          continue;
//
//        for (int j = 0; j < this->phi->phi_supp[i].size(); ++j) {
//          if (this->phi->phi_supp[i].at(j) == (this->M - 1)) {
//            this->phi->phi_supp[i].set(j, idx);
//          }
//        }
//        for (int j = 0; j < this->phi->phi_inv_supp[i].size(); ++j) {
//          if (this->phi->phi_inv_supp[i].at(j) == (this->M - 1)) {
//            this->phi->phi_inv_supp[i].set(j, idx);
//          }
//        }
//      }
//
//      this->phi->phi_supp[hap_id] = this->phi->phi_supp[this->M - 1];
//      this->phi->phi_inv_supp[hap_id] = this->phi->phi_inv_supp[this->M - 1];
//      this->phi->phi_supp.pop_back();
//      this->phi->phi_inv_supp.pop_back();
//      this->phi->phi_vec[hap_id] = this->phi->phi_vec[this->M - 1];
//      this->phi->phi_inv_vec[hap_id] = this->phi->phi_inv_vec[this->M - 1];
//      this->phi->phi_vec.pop_back();
//      this->phi->phi_inv_vec.pop_back();
//
//      // Update pref
//      for (int col = 0; col < this->N; ++col) {
//        // Update pref beg
//        for (int i = 0; i < this->columns[col].pref_samples_beg.size(); ++i) {
//          if (columns[col].pref_samples_beg.at(i) == (this->M - 1)) {
//            columns[col].pref_samples_beg.set(i, hap_id);
//          }
//        }
//        // Update pref end
//        for (int i = 0; i < this->columns[col].pref_samples_end.size(); ++i) {
//          if (columns[col].pref_samples_end.at(i) == (this->M - 1)) {
//            columns[col].pref_samples_end.set(i, hap_id);
//          }
//        }
//      }
//    } else {
//      this->phi->phi_inv_vec.pop_back();
//      this->phi->phi_vec.pop_back();
//      this->phi->phi_supp.pop_back();
//      this->phi->phi_inv_supp.pop_back();
//    }

    haplotype_ids.erase(hap_id);
    --this->M;
    --this->phi->total_haplotypes;
  }

  // Build the reference panel
  void BuildFromVCF(std::string &filename, bool verbose) {
    std::string line = "##";
    if (std::ifstream inFile(filename); inFile.is_open()) {
      // go through headers
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
      vector<unsigned int> prefix_arr(M, 0);
      std::iota(prefix_arr.begin(), prefix_arr.end(), 0);
      // TODO: Possible optimization using sdsl int vec
      vector<vector<unsigned int>> sites_where_sample_beg(M);
      vector<vector<unsigned int>> sites_where_sample_end(M);
      vector<unsigned int> div(M, 0);
      vector<unsigned int> div_u, div_v;
      unsigned int p = col + 1, q = col + 1;

      while (getline(inFile, line)) {
        std::istringstream iss(line);
        token = "";

        // read through the line and store in single_col
        std::vector<bool> single_col;
        for (int i = 0; i < (this->M / 2) + 9; ++i) {
          iss >> token;
          if (i < 9) {
            continue;
          }
          single_col.push_back(static_cast<bool>(token[0] - '0'));
          single_col.push_back(static_cast<bool>(token[2] - '0'));
        }

        packed_spsi temp_zeros;
        packed_spsi temp_ones;
        packed_spsi temp_combined;
//        packed_spsi temp_sample_beg;
//        packed_spsi temp_sample_end;
        succinct_spsi temp_sample_beg;
        succinct_spsi temp_sample_end;
        succinct_spsi temp_div;
        bool start_with_zero = false;

        assert(single_col.size() == this->M);
        p = col + 1;
        q = col + 1;
        for (int i = 0; i < this->M; ++i) {
          // first allele
          if (i == 0) {
            if (single_col[prefix_arr[i]]) { // allele: 1
              v.push_back(prefix_arr[i]);
              div_v.push_back(q);
              q = 0;
            } else { // allele: 0
              u.push_back(prefix_arr[i]);
              start_with_zero = true;
              div_u.push_back(p);
              p = 0;
            }
            temp_sample_beg.push_back(prefix_arr[i]);
            temp_div.push_back(div[i]);
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
            temp_div.push_back(div[i]);
          } else {
            ++cnt;
          }

          if (div[i] > p) {
            p = div[i];
          }
          if (div[i] > q) {
            q = div[i];
          }

          if (single_col[prefix_arr[i]]) { // allele 1
            v.push_back(prefix_arr[i]);
            div_v.push_back(q);
            q = 0;
          } else { // allele 0
            u.push_back(prefix_arr[i]);
            div_u.push_back(p);
            p = 0;
          }
        }
        temp_sample_end.push_back(prefix_arr[M - 1]);


        // edge case
        if (cnt > 0)
          freq.push_back(cnt);

        // populate dynamic data structures
        for (int i = 0; i < freq.size(); ++i) {
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
        assert(temp_div.size() == temp_combined.size());

        for (auto i = 0; i < temp_sample_beg.size(); ++i) {
          sites_where_sample_beg[temp_sample_beg.at(i)].push_back(col);
          sites_where_sample_end[temp_sample_end.at(i)].push_back(col);
        }

        // build each column
        dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                           (std::move(temp_sample_beg)), (std::move(temp_sample_end)), (std::move(temp_div)),
                           start_with_zero, u.size());
        columns.emplace_back(coln);

        total_runs += freq.size();

        // next col prefix arr
        prefix_arr.clear();
        prefix_arr.insert(prefix_arr.end(), u.begin(), u.end());
        prefix_arr.insert(prefix_arr.end(), v.begin(), v.end());
        div.clear();
        div.insert(div.end(), div_u.begin(), div_u.end());
        div.insert(div.end(), div_v.begin(), div_v.end());
        div_u.clear();
        div_v.clear();
        ++col;
        u.clear();
        v.clear();
        freq.clear();
      }
      assert(columns.size() == N);
      // build phi data-structure
      this->phi = new phi_ds(columns, M, N, sites_where_sample_beg, sites_where_sample_end, prefix_arr, div, verbose);
      assert(col == N);
      total_runs += freq.size();
      cout << "Phi support size (in bytes) = " << this->phi->size_in_bytes(verbose) << "\n";
      cout << "Avg runs = " << static_cast<float>(total_runs) / N << "\n";

      inFile.close();
    } else {
      std::cerr << "Couldn't find : " << filename << "\n";
      exit(1);
    }

  }

  /*
  void Build(std::vector<std::vector<bool>> &alleles) {

    vector<int> u, v;
    vector<int> freq;
    int total_runs = 0;
    int col = 0;
    int cnt = 1;
    bool prev_allele = false;
    vector<int> prefix_arr(M, 0);
    std::iota(prefix_arr.begin(), prefix_arr.end(), 0);

    vector<vector<unsigned int>> sites_where_sample_beg(M);
    vector<vector<unsigned int>> sites_where_sample_end(M);

    while (col < N) {
      packed_spsi temp_zeros;
      packed_spsi temp_ones;
      packed_spsi temp_combined;
//      packed_spsi temp_sample_beg;
//      packed_spsi temp_sample_end;
      succinct_spsi temp_sample_beg;
      succinct_spsi temp_sample_end;
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
      temp_sample_end.push_back(prefix_arr[M - 1]);

      // edge case
      if (cnt > 0)
        freq.push_back(cnt);

      // populate dynamic data structures
      for (int i = 0; i < freq.size(); ++i) {
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

      for (auto i = 0; i < temp_sample_beg.size(); ++i) {
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
    cout << "Avg runs = " << static_cast<float>(total_runs) / N << "\n";
  }
   */
};
