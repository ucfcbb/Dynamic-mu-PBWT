#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <sstream>
#include <exception>

#include "dynamic/dynamic.hpp"
#include "phi.h"
#include "dcpbwt_column.h"
#include "utils.h"

#include "htslib/vcf.h"

using namespace dyn;
using namespace std;

class DCPBWT {
    public:
    unsigned int M = 0; // #haplotypes
    unsigned int N = 0; // #sites
    std::vector<dcpbwt_column> columns;
    phi_ds phi;

    // default constructor
    DCPBWT() : M(0), N(0) {}

    // constructor
    DCPBWT(std::string ref_vcf_file) {
//     Build(ref_vcf_file.c_str(), verbose);
        BuildFromVCF(ref_vcf_file);
        std::cout << "Total haplotypes: " << this->M << std::endl;
        std::cout << "Total sites: " << this->N << std::endl;
    }
    ~DCPBWT() {
    }


    // Helper Methods
    /*
     * Returns the index of run that the haplotype index i belongs to.
     */
    [[nodiscard]] unsigned int get_run_idx(const unsigned int col, const unsigned int i) const {
        if (i >= (this->M - 1)) {
            return this->columns[col].combined.size() - 1;
        }
        return this->columns[col].combined.search(i + 1);
    }

    /*
     * Returns the index of the head of a run specified by run_idx
     */
    [[nodiscard]] unsigned int get_run_head(const unsigned int col, const unsigned int run_idx) const {
        if (run_idx >= this->columns[col].combined.size()) {
            cerr << "Out of bounds: Accessing run that doesn't exist!\n";
            exit(EXIT_FAILURE);
        }
        if (run_idx == 0) {
            return 0;
        }
        return this->columns[col].combined.psum(run_idx - 1);
    }

    /*
     * Returns the allele value of the run specified by run_idx
     */
    [[nodiscard]] bool get_run_val(const unsigned int col, const unsigned int run_idx) const {
        if (!(run_idx & 1)) {
            if (this->columns[col].start_with_zero) {
                return false;
            }
            return true;
        }
        if (this->columns[col].start_with_zero) {
            return true;
        }
        return false;
    }

    /*
     * Returns true if at a run end
     */
    [[nodiscard]] bool isRunEnd(const unsigned int col, const unsigned int run_idx, const unsigned int hap_id) const {
        return hap_id == this->columns[col].pref_samples_end.at(run_idx);
    }

    /*
     * Returns true if at a run start
     */
    [[nodiscard]] bool isRunStart(const unsigned int col, const unsigned int run_idx, const unsigned int hap_id) const {
        return hap_id == this->columns[col].pref_samples_beg.at(run_idx);
    }

    /*
     * Returns true if a run has a single element
     */
    [[nodiscard]] bool isRunWithSingleElement(const unsigned int col, const unsigned int run_idx) const {
        return this->columns[col].combined.at(run_idx) == 1;
    }

    /*
     * Maps the haplotype at the bottom of a column k to column k+1
     */
    [[nodiscard]] std::pair<unsigned int, unsigned int> column_end_lf(unsigned int col, bool val) const {
        if (val) { return {this->M, UINT_MAX}; } // last row, invalid hapID
        // Only zeros no ones in a column
        if (!val && this->columns[col].ones.size() == 0) { return {this->M, UINT_MAX}; } // last row, invalid hapID
        unsigned int ind = this->columns[col].num_zeros;
        unsigned int ridx = 0; // assume column starts with a run of ones
        if (this->columns[col].start_with_zero) { // the second run at index 1 will be run of ones
            ridx = 1;
        }
        return {ind, this->columns[col].pref_samples_beg.at(ridx)};
    }

    /*
     * modified extension that maps i-th haplotype in col to col+1
     */
    std::pair<unsigned int, unsigned int> w_mod(const unsigned int i,
                                                const unsigned int col,
                                                const bool val,
                                                const unsigned int pref_val) {
        if (i == this->M) { return column_end_lf(col, val); }
        const unsigned int run_idx = get_run_idx(col, i);
        const bool rval = get_run_val(col, run_idx);
        if (rval == val) {
            return {lf(col, i, val, run_idx, rval), pref_val};
        }

        if (run_idx + 1 >= this->columns[col].combined.size())
            return column_end_lf(col, val);
        return {lf(col, get_run_head(col, run_idx + 1), val, run_idx + 1, get_run_val(col, run_idx + 1)),
                this->columns[col].pref_samples_beg[run_idx + 1]}; // sample_beg fetch pref-val
    }

    /*
     * FL mapping subroutine for w_mod
     */
    unsigned int lf(const unsigned int col,
                    const unsigned int i,
                    const bool value,
                    unsigned int run_idx,
                    unsigned int rval) {
        if (i == M) {
            if (value) // val is 1
                return M;
            return columns[col].num_zeros; // val is 0
        }

        //TODO: this is redundant and should be passed as an argument
        // unsigned int run_idx = get_run_idx(col, i);
        unsigned int run_head = get_run_head(col, run_idx);
        unsigned int offset = i - run_head;

        if (value != rval) {
//        if (value != get_run_val(col, run_idx)) {
            offset = 0;
        }

        auto uv = uv_trick(col, run_idx);
//      auto uv = uv_trick(col, i);

        if (value) {
            return columns[col].num_zeros + uv.second + offset;
        }
        return uv.first + offset;
    }

    /*
     * Reverse FL mapping from col k+1 to col k
     */
    unsigned int reverse_lf(const unsigned int col, unsigned int idx) {
        if (col == 0) {
            return 0;
        }

        unsigned int prev_col = col - 1;
        unsigned int num_zeros = this->columns[prev_col].num_zeros;
        // came from a run of zeros in prev col
        unsigned offset = 0;
        if (idx < num_zeros) {
            unsigned int zero_idx = this->columns[prev_col].zeros.search(idx + 1);
            unsigned int zeros_before = 0;
            if (zero_idx > 0)
                zeros_before = this->columns[prev_col].zeros.psum(zero_idx - 1);
            offset = idx - zeros_before;

            // find run index
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
        } else { // came from a run of ones
            idx -= num_zeros;

            unsigned int ones_idx = this->columns[prev_col].ones.search(idx + 1);

            // calculate run idx
            unsigned int run_idx = 0;
            if (this->columns[prev_col].start_with_zero) {
                run_idx = ones_idx * 2 + 1;
            } else {
                run_idx = ones_idx * 2;
            }

            // calculate offset
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

    /*
     * Returns the number of zeros before a given run_idx
     */
    unsigned int zeros_before(const unsigned int col, const unsigned int run_idx) {
        unsigned int retval = 0;
        assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
        if (columns[col].start_with_zero) {
            if (run_idx % 2 == 0)
                retval = this->columns[col].zeros.psum((run_idx - 2) / 2);
            else
                retval = this->columns[col].zeros.psum((run_idx - 1) / 2);
        } else {
            if (run_idx % 2 == 0)
                retval = this->columns[col].zeros.psum((run_idx - 1) / 2);
            else
                retval = this->columns[col].zeros.psum((run_idx - 2) / 2);
        }
        return retval;
    }

    /*
     * Returns the number of ones before a given run_idx
     */
    unsigned int ones_before(const unsigned int col, const unsigned int run_idx) {
        unsigned int retval = 0;
        assert(((run_idx - 1) >= 0) && ((run_idx - 1) < UINT_MAX));
        if (columns[col].start_with_zero) {
            if (run_idx % 2 == 0)
                retval = this->columns[col].ones.psum((run_idx - 1) / 2);
            else
                retval = this->columns[col].ones.psum((run_idx - 2) / 2);
        } else {
            if (run_idx % 2 == 0)
                retval = this->columns[col].ones.psum((run_idx - 2) / 2);
            else
                retval = this->columns[col].ones.psum((run_idx - 1) / 2);
        }
        return retval;
    }


    /*
     * Returns the number of zeros and ones before ith haplotype in a column
     */
//  pair<unsigned int, unsigned int> uv_trick(const unsigned int col, const unsigned i) {
    // unsigned int run_idx = get_run_idx(col, i);
    pair<unsigned int, unsigned int> uv_trick(const unsigned int col, const unsigned run_idx) {
        unsigned int u = 0, v = 0;
        if (run_idx == 0) {
            return make_pair(0, 0);
        }
        if (run_idx == 1) {
            if (this->columns[col].start_with_zero) {
                return make_pair(zeros_before(col, run_idx), 0);
            }
            return make_pair(0, ones_before(col, run_idx));
        }
        u = zeros_before(col, run_idx);
        v = ones_before(col, run_idx);
        return make_pair(u, v);
    }

    /*
     * Returns true if allele matches the value of run indexed by run_idx
     */
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

    /*
     * Inserting at the start of a run
     */
    void InsertAtRunStart(const unsigned int col,
                          const unsigned idx,
                          const unsigned int hap_id,
                          const unsigned int inserted_hap_id,
                          bool allele,
                          packed_spsi &temp_supp,
                          packed_spsi &temp_inv_supp,
                          packed_spsi &temp_div_supp,
                          vector<unsigned int> &temp_div_query,
                          vector<unsigned int> &temp_div_below_query) {

        unsigned int run_idx = get_run_idx(col, idx);
        unsigned int hap_before = 0;
        if (idx == 0) {
            hap_before = inserted_hap_id;
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
            this->columns[col].pref_samples_beg.set(run_idx, inserted_hap_id);

            // unset old head of run
            auto col_rank = this->phi.phi_vec[hap_id].rank1(col);
            this->phi.phi_vec[hap_id].set(col, false);
            // Handle for the N-th column
            if (col == this->N) {
                this->phi.phi_supp[hap_id].set(this->phi.phi_supp[hap_id].size() - 1, inserted_hap_id);
                this->phi.phi_supp_lcp[hap_id].set(this->phi.phi_supp_lcp[hap_id].size() - 1,
                                                   temp_div_below_query[col]);
                temp_inv_supp.push_back(hap_id);
            } else {
                this->phi.phi_supp[hap_id].remove(col_rank);
                this->phi.phi_supp_lcp[hap_id].remove(col_rank);
            }

            // set phi for new run head
            this->phi.phi_vec[inserted_hap_id].set(col, true);
            temp_supp.push_back(hap_before);

            // update phi for hap_before as well
            if (hap_before != inserted_hap_id) {
                col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
                this->phi.phi_inv_supp[hap_before].set(col_rank, inserted_hap_id);
            }

            // Update div sample
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
                this->columns[col].pref_samples_beg.insert(run_idx, inserted_hap_id);
                this->columns[col].pref_samples_end.insert(run_idx, inserted_hap_id);

                // update phi for old run head
                auto col_rank = this->phi.phi_vec[hap_id].rank1(col);
                this->phi.phi_supp[hap_id].set(col_rank, inserted_hap_id);
                this->phi.phi_supp_lcp[hap_id].set(col_rank, temp_div_below_query[col]);

                // update for newly inserted head
                this->phi.phi_vec[inserted_hap_id].set(col, true);
                this->phi.phi_inv_vec[inserted_hap_id].set(col, true);
                temp_supp.push_back(inserted_hap_id);
                temp_inv_supp.push_back(hap_id);

                // Update div sample
                // Insert new div sample and also update the one that's below it
                // Update the one that's below first
                this->columns[col].div_samples_beg.set(run_idx, temp_div_below_query[col]);
                assert(temp_div_query[col] == col);
                this->columns[col].div_samples_beg.insert(0, temp_div_query[col]);
                temp_div_supp.push_back(temp_div_query[col]);
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
                this->columns[col].pref_samples_end.set(run_idx - 1, inserted_hap_id);

                // update phi for old hap at prev run's end
                auto col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
                if (col == this->N) {
                    this->phi.phi_inv_supp[hap_before].set(col_rank, inserted_hap_id);
                    temp_supp.push_back(hap_before);
                    temp_div_supp.push_back(temp_div_query[col]);
                } else {
                    this->phi.phi_inv_supp[hap_before].remove(col_rank);
                }
                this->phi.phi_inv_vec[hap_before].set(col, false);

                // update for newly inserted hap
                this->phi.phi_inv_vec[inserted_hap_id].set(col, true);
                temp_inv_supp.push_back(hap_id);

                // update phi_inv value for hap_id
                // if (hap_id != UINT_MAX)  <- this should NEVER happen
                assert(hap_id != inserted_hap_id);
                col_rank = this->phi.phi_vec[hap_id].rank1(col);
                this->phi.phi_supp[hap_id].set(col_rank, inserted_hap_id);
                // Update div sample of the seq below
                this->phi.phi_supp_lcp[hap_id].set(col_rank, temp_div_below_query[col]);
                this->columns[col].div_samples_beg.set(run_idx, temp_div_below_query[col]);
            }
        }
    }

    /*
     * Inserting at the bottom of a column i.e. the last run in the column
     */
    void InsertAtBottom(const unsigned int col,
                        const unsigned int run_idx,
                        const unsigned int hap_id,
                        const unsigned int inserted_hap_id,
                        const bool allele,
                        packed_spsi &temp_supp,
                        packed_spsi &temp_inv_supp,
                        packed_spsi &temp_div_supp,
                        vector<unsigned int> &temp_div_query) {
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
            this->columns[col].pref_samples_end.set(run_idx, inserted_hap_id);

            // Update old bottom's Phi_inv structures
            auto col_rank = this->phi.phi_inv_vec[old_bottom].rank1(col);
            if (col == this->N) {
                this->phi.phi_inv_supp[old_bottom].set(col_rank, inserted_hap_id);
                temp_supp.push_back(old_bottom);
                temp_div_supp.push_back(temp_div_query[col]);
            } else {
                this->phi.phi_inv_supp[old_bottom].remove(col_rank);
            }
            this->phi.phi_inv_vec[old_bottom].set(col, false);

            // Update new bottom's Phi_inv structures
            this->phi.phi_inv_vec[inserted_hap_id].set(col, true);
            temp_inv_supp.push_back(inserted_hap_id);
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

            // Update pref beg/end
            this->columns[col].pref_samples_beg.push_back(inserted_hap_id);
            this->columns[col].pref_samples_end.push_back(inserted_hap_id);

            // Update phi supports for inserted haplotype
            this->phi.phi_vec[inserted_hap_id].set(col, true);
            this->phi.phi_inv_vec[inserted_hap_id].set(col, true);
            temp_supp.push_back(hap_before);
            temp_inv_supp.push_back(inserted_hap_id);

            // Update phi_inv for hap_before
            auto col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
            this->phi.phi_inv_supp[hap_before].set(col_rank, inserted_hap_id);

            // Insert Div val since this is a new run
            this->columns[col].div_samples_beg.push_back(temp_div_query[col]);
            temp_div_supp.push_back(temp_div_query[col]);
        }
    }

    /*
     * Insert in a single column:
     *  hap_id: hap_id of the sequence above which we insert inserted_hap_id
     *  idx: index of hap_id in col
     *  allele: allele of the inserted_hap_id bit at col
     */
    void Insert(const unsigned int col,
                const unsigned idx,
                const unsigned int hap_id,
                const unsigned int inserted_hap_id,
                bool allele,
                packed_spsi &temp_supp,
                packed_spsi &temp_inv_supp,
                packed_spsi &temp_div_supp,
                vector<unsigned int> &temp_div_query,
                vector<unsigned int> &temp_div_below_query) {
        // Handle when panel is empty
        if (this->M == 0) {
            assert(this->columns[col].combined.size() == 0);
            assert(this->columns[col].ones.size() == 0);
            assert(this->columns[col].zeros.size() == 0);
            assert(this->columns[col].num_zeros == 0);
            this->columns[col].combined.push_back(1);
            if (allele) {
                this->columns[col].ones.push_back(1);
                this->columns[col].start_with_zero = false;
            } else {
                this->columns[col].zeros.push_back(1);
                ++this->columns[col].num_zeros;
                this->columns[col].start_with_zero = true;
            }
            // Update pref beg/end
            assert(this->columns[col].pref_samples_beg.size() == 0);
            this->columns[col].pref_samples_beg.push_back(inserted_hap_id);
            assert(this->columns[col].pref_samples_end.size() == 0);
            this->columns[col].pref_samples_end.push_back(inserted_hap_id);
            assert(this->columns[col].div_samples_beg.size() == 0);
            this->columns[col].div_samples_beg.push_back(col);

            // Update div beg
            temp_div_supp.push_back(col);
            temp_supp.push_back(inserted_hap_id);
            temp_inv_supp.push_back(inserted_hap_id);
            return;
        }

        auto run_idx = get_run_idx(col, idx);
        if (isRunStart(col, run_idx, hap_id)) {
            InsertAtRunStart(col,
                             idx,
                             hap_id,
                             inserted_hap_id,
                             allele,
                             temp_supp,
                             temp_inv_supp,
                             temp_div_supp,
                             temp_div_query,
                             temp_div_below_query);
            return;
        }

        if (idx == this->M) {
            InsertAtBottom(col,
                           run_idx,
                           hap_id,
                           inserted_hap_id,
                           allele,
                           temp_supp,
                           temp_inv_supp,
                           temp_div_supp,
                           temp_div_query);
            return;
        }

        // if run value is same as the allele
        if (AlleleMatchesRun(col, run_idx, allele)) {
            this->columns[col].combined.increment(run_idx, 1);
            if (allele) {
                this->columns[col].ones.increment(run_idx / 2, 1);
            } else {
                this->columns[col].zeros.increment(run_idx / 2, 1);
                ++columns[col].num_zeros;
            }

            if (col == this->N) {
                auto hap_before_opt = this->phi.phi(hap_id, col);
                assert(hap_before_opt.has_value());
                unsigned int hap_before = hap_before_opt.value();
                unsigned int hap_after = hap_id;

                this->phi.phi_supp[hap_after].set(this->phi.phi_supp[hap_after].size() - 1, inserted_hap_id);
                this->phi.phi_supp_lcp[hap_after].set(this->phi.phi_supp_lcp[hap_after].size() - 1,
                                                      temp_div_below_query[col]);
                this->phi.phi_inv_supp[hap_before].set(this->phi.phi_inv_supp[hap_before].size() - 1, inserted_hap_id);
                temp_supp.push_back(hap_before);
                temp_inv_supp.push_back(hap_after);
                temp_div_supp.push_back(temp_div_query[col]);
            }
        } else {
            /*
             * Create a new run as the existing run will be split
             */
            auto hap_before_opt = this->phi.phi(hap_id, col);
            assert(hap_before_opt.has_value());
            unsigned int hap_before = hap_before_opt.value();
            unsigned int hap_after = hap_id;

            // Update run info
            unsigned int run_head = get_run_head(col, run_idx);
            unsigned int left_rem = idx - run_head;
            unsigned int right_rem = columns[col].combined.at(run_idx) - left_rem;
            this->columns[col].combined.set(run_idx, left_rem);
            this->columns[col].combined.insert(run_idx + 1, 1);
            this->columns[col].combined.insert(run_idx + 2, right_rem);

            if (allele) {
                // Inserting 1 in a run of 0s
                if (run_idx == 0) {
                    this->columns[col].ones.insert(run_idx, 1);
                } else {
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
                this->columns[col].pref_samples_beg.push_back(inserted_hap_id);
                this->columns[col].pref_samples_beg.push_back(hap_after);
            } else {
                this->columns[col].pref_samples_beg.insert(run_idx + 1, inserted_hap_id);
                this->columns[col].pref_samples_beg.insert(run_idx + 2, hap_after);
            }
            // Update Pref end
            this->columns[col].pref_samples_end.insert(run_idx, inserted_hap_id);
            this->columns[col].pref_samples_end.insert(run_idx, hap_before);

            // Update phi supports for hap_before
            auto col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
            this->phi.phi_inv_supp[hap_before].insert(col_rank, inserted_hap_id);
            this->phi.phi_inv_vec[hap_before].set(col, true);

            // Update phi supports for hap_after
            col_rank = this->phi.phi_vec[hap_after].rank1(col);
            this->phi.phi_supp[hap_after].insert(col_rank, inserted_hap_id);
            this->phi.phi_vec[hap_after].set(col, true);
            this->phi.phi_supp_lcp[hap_after].insert(col_rank, temp_div_below_query[col]);

            // update phi structure for Inserted Haplotype
            this->phi.phi_vec[inserted_hap_id].set(col, true);
            this->phi.phi_inv_vec[inserted_hap_id].set(col, true);
            temp_supp.push_back(hap_before);
            temp_inv_supp.push_back(hap_after);

            // Insert div value for inserted hap and the hap below it since we're forming two new head of runs
            if (run_idx + 1 >= this->columns[col].combined.size()) {
                this->columns[col].div_samples_beg.push_back(temp_div_query[col]);
                this->columns[col].div_samples_beg.push_back(temp_div_below_query[col]);

                this->phi.phi_supp_lcp[hap_after].push_back(temp_div_below_query[col]);
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

    void CalculateDivergenceVal(vector<pair<unsigned int, unsigned int>> &insertion_indices,
                                vector<bool> &insertion_alleles,
                                vector<unsigned int> &temp_div_seqn,
                                vector<unsigned int> &temp_div_below_seqn) {
        unsigned int zs = 0; // starting position of match between inserted seqn and seqn above it
        unsigned int bs = 0; // starting position of match between inserted seqn and seqn below it
        unsigned int hap_id_above_prev = -1;
        for (unsigned int k = this->N; k > 0; --k) {
            zs = k;
            bs = k;
            unsigned int position = insertion_indices[k].first;
            unsigned int hap_id = insertion_indices[k].second;

            // update divergence values for sequence BELOW query
            // "starting pos" definition of divergence value
            if (k != this->N && insertion_indices[k + 1].second == hap_id && temp_div_below_seqn[k + 1] < k)
                bs = temp_div_below_seqn[k + 1];
            else if (position != this->M) {
                while (bs > 0 &&
                    get_curr_char(position, bs) == insertion_alleles[bs - 1]) {
                    position = reverse_lf(bs, position);
                    --bs;
                }
            }

            // Update divergence values for query sequence
            // "starting pos" definition of divergence value
            position = insertion_indices[k].first;
            if (position-- != 0) {
                unsigned int hap_id_above;// = this->phi->phi(hap_id, k);;
                if (position == this->M - 1) {
                    hap_id_above = this->columns[k].pref_samples_end.at(this->columns[k].combined.size() - 1);
                } else {
                    auto hap_id_above_opt = this->phi.phi(hap_id, k);
                    assert(hap_id_above_opt.has_value());
                    hap_id_above = hap_id_above_opt.value();
                }
                if (k != this->N && hap_id_above_prev == hap_id_above && temp_div_seqn[k + 1] < k)
                    zs = temp_div_seqn[k + 1];
                else {
                    while (zs > 0 &&
                        get_curr_char(position, zs) == insertion_alleles[zs - 1]) {
                        position = reverse_lf(zs, position);
                        --zs;
                    }
                }
                hap_id_above_prev = hap_id_above;
            } else {
                hap_id_above_prev = this->M;
            }
            temp_div_seqn[k] = zs;
            temp_div_below_seqn[k] = bs;
        }
        temp_div_seqn[0] = temp_div_below_seqn[0] = 0;
    }

    void InsertFirstHaplotype(std::vector<bool> &query) {
        unsigned int cnt_0 = 0;
        for (unsigned int i = 0; i < query.size(); ++i) {
            packed_spsi temp_ones, temp_zeros, temp_combined;
            packed_spsi temp_sample_beg, temp_sample_end, temp_div;
            cnt_0 = 0;
            if (query[i]) {
                temp_ones.push_back(1);
            } else {
                temp_zeros.push_back(1);
                cnt_0 = 1;
            }
            temp_combined.push_back(1);
            temp_sample_beg.push_back(this->M);
            temp_sample_end.push_back(this->M);
            temp_div.push_back(i);
            // build each column
            dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                               (std::move(temp_sample_beg)), (std::move(temp_sample_end)), (std::move(temp_div)),
                               !query[i], cnt_0);
            this->columns.push_back(coln);
        }
        // Handle for N-th coln
        packed_spsi temp_ones, temp_zeros, temp_combined;
        packed_spsi temp_sample_beg, temp_sample_end, temp_div;
        if (query[query.size() - 1]) {
            temp_ones.push_back(1);
        } else {
            temp_zeros.push_back(1);
            cnt_0 = 1;
        }
        temp_combined.push_back(1);
        temp_sample_beg.push_back(this->M);
        temp_sample_end.push_back(this->M);
        temp_div.push_back(query.size());
        dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                           (std::move(temp_sample_beg)), (std::move(temp_sample_end)), (std::move(temp_div)),
                           !query[query.size() - 1], cnt_0);
        this->columns.push_back(coln);
        assert(this->columns.size() == query.size() + 1);

        // initialize new phi structure
        this->phi = phi_ds();
        suc_bv tmp_b, tmp_e;
        packed_spsi temp_supp;
        packed_spsi temp_inv_supp;
        packed_spsi temp_div_supp;
        for (unsigned int col = 0; col <= query.size(); ++col) {
            tmp_b.push_back(true);
            tmp_e.push_back(true);
            temp_supp.push_back(this->M);
            temp_inv_supp.push_back(this->M);
            temp_div_supp.push_back(col);
        }
        this->phi.phi_vec.push_back(tmp_b);
        this->phi.phi_inv_vec.push_back(tmp_e);
        this->phi.phi_supp.push_back(temp_supp);
        this->phi.phi_inv_supp.push_back(temp_inv_supp);
        this->phi.phi_supp_lcp.push_back(temp_div_supp);
        ++this->M;
        ++this->phi.total_haplotypes;
        this->N = query.size();
    }

    void InsertSingleHaplotype(std::vector<bool> &query) {
        if (this->M == 0) {
            InsertFirstHaplotype(query);
            return;
        }

        assert(query.size() == N);
        vector<pair<unsigned int, unsigned int>> insertion_indices(this->N + 1); // stores {index, hapid}
        insertion_indices[0].first = this->M;
        insertion_indices[0].second = this->M;

        // Calculate insertion indices
        for (unsigned int col = 0; col < query.size(); ++col) {
            insertion_indices[col + 1] =
                w_mod(insertion_indices[col].first, col, query[col], insertion_indices[col].second);
        }

        // Calculate divergence values
        vector<unsigned int> temp_div_seqn(this->N + 1, 0);
        vector<unsigned int> temp_div_below_seqn(this->N + 1, 0);
        CalculateDivergenceVal(insertion_indices, query, temp_div_seqn, temp_div_below_seqn);

        // initialize new phi structure
        suc_bv tmp_b, tmp_e;
        for (unsigned int col = 0; col <= query.size(); ++col) {
            tmp_b.push_back(false);
            tmp_e.push_back(false);
        }
        this->phi.phi_vec.push_back(tmp_b);
        this->phi.phi_inv_vec.push_back(tmp_e);
        packed_spsi temp_supp;
        packed_spsi temp_inv_supp;
        packed_spsi temp_div_supp;

        // perform insertion
        // where div samples are also updated accordingly
        for (unsigned int col = 0; col <= query.size(); ++col) {
            if (col == query.size()) {
                Insert(col,
                       insertion_indices[col].first,
                       insertion_indices[col].second,
                       this->M,
                       query[col - 1],
                       temp_supp,
                       temp_inv_supp,
                       temp_div_supp,
                       temp_div_seqn,
                       temp_div_below_seqn);

            } else {
                Insert(col,
                       insertion_indices[col].first,
                       insertion_indices[col].second,
                       this->M,
                       query[col],
                       temp_supp,
                       temp_inv_supp,
                       temp_div_supp,
                       temp_div_seqn,
                       temp_div_below_seqn);
            }
        }

        assert(temp_supp.size() == temp_div_supp.size());
        // handling for col == N case
        this->phi.phi_supp.push_back(temp_supp);
        this->phi.phi_inv_supp.push_back(temp_inv_supp);
        this->phi.phi_supp_lcp.push_back(temp_div_supp);
        // increment total # of haplotypes
        ++this->M;
        ++this->phi.total_haplotypes;
    }

    /*
     * Delete run with a single value remaining. The adjacent runs (if exist) will be merged in this case.
     */
    void DeleteSingleRun(const unsigned int col, const unsigned int idx, const unsigned int hap_id, const bool allele) {
        unsigned int run_idx = 0;

        // Handle deletion of the first run
        if (idx == 0) {
            // new top will be the head of next run
            unsigned int new_top = hap_id;

            // check if it's not the last run i.e. single haplotype remaining
            if (idx + 1 < this->columns[col].pref_samples_beg.size()) {
                assert(this->M > 1);
                new_top = this->columns[col].pref_samples_beg.at(idx + 1);
            }

            // Update pref_beg/end
            this->columns[col].pref_samples_beg.remove(run_idx);
            this->columns[col].pref_samples_end.remove(run_idx);

            // Update run info
            this->columns[col].combined.remove(idx);
            if (allele) {
                this->columns[col].ones.remove(idx);
                this->columns[col].start_with_zero = true;
            } else {
                this->columns[col].zeros.remove(idx);
                this->columns[col].start_with_zero = false;
                --this->columns[col].num_zeros;
            }

            // TODO: This clean up can be performed at the very end and will probably be easier
            // Update phi structures for removed haplotype
            auto col_rank = this->phi.phi_vec[hap_id].rank1(col);
            // so that the phi_supp info in the N-th column isn't deleted
            if (this->phi.phi_supp[hap_id].size() > 1) {
                this->phi.phi_supp[hap_id].remove(col_rank);
                this->phi.phi_supp_lcp[hap_id].remove(col_rank);
            }
            // Update phi_inv
            col_rank = this->phi.phi_inv_vec[hap_id].rank1(col);
            // so that the phi_inv_supp info in the N-th column isn't deleted
            if (this->phi.phi_inv_supp[hap_id].size() > 1)
                this->phi.phi_inv_supp[hap_id].remove(col_rank);
            this->phi.phi_vec[hap_id].set(col, false);
            this->phi.phi_inv_vec[hap_id].set(col, false);

            // Remove div sample of the removed haplotype
            this->columns[col].div_samples_beg.remove(idx);

            // If not the last haplotype, set the new top's div/phi value
            if (new_top != hap_id) {
                assert(this->M > 1);
                this->columns[col].div_samples_beg.set(idx, col);
                // update phi for the new top
                col_rank = this->phi.phi_vec[new_top].rank1(col);
                this->phi.phi_supp[new_top].set(col_rank, new_top);
                this->phi.phi_supp_lcp[new_top].set(col_rank, col);
            }
            return;
        }

        // Handle deletion of the bottom run in a column
        if (idx == this->M - 1) {
            assert(this->M > 1);
            run_idx = this->columns[col].pref_samples_beg.size() - 1;

            // new bottom will be the bottom of previous run
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

            // TODO: This can be taken care of at the very end and will be easier to clean up
            // update phi structures for removed haplotype
            auto col_rank = this->phi.phi_vec[hap_id].rank1(col);
            if (this->phi.phi_supp[hap_id].size() > 1) {
                this->phi.phi_supp[hap_id].remove(col_rank);
                this->phi.phi_supp_lcp[hap_id].remove(col_rank);
            }
            col_rank = this->phi.phi_inv_vec[hap_id].rank1(col);
            if (this->phi.phi_inv_supp[hap_id].size() > 1)
                this->phi.phi_inv_supp[hap_id].remove(col_rank);
            this->phi.phi_vec[hap_id].set(col, false);
            this->phi.phi_inv_vec[hap_id].set(col, false);
            // Update div sample
            this->columns[col].div_samples_beg.remove(run_idx);

            // Update phi_inv for the new bottom
            col_rank = this->phi.phi_inv_vec[new_bottom].rank1(col);
            this->phi.phi_inv_supp[new_bottom].set(col_rank, new_bottom);

            return;
        }

        // Delete run from the middle of the column
        run_idx = get_run_idx(col, idx);
        unsigned int hap_before = this->columns[col].pref_samples_end[run_idx - 1];
        unsigned int hap_after = this->columns[col].pref_samples_beg[run_idx + 1];
        //unsigned int hap_id_div_val = this->columns[col].div_samples_beg.at(run_idx);
        //unsigned int hap_after_div_val = this->columns[col].div_samples_beg.at(run_idx + 1);

        /* 0 0 0 0 0 0 0  |  1   | 0 0 0 0 0 0
         *    ridx-1      | ridx | ridx + 1
         * merge ridx + 1 into ridx - 1
         */

        // Update run info spsis
        this->columns[col].combined.increment(run_idx - 1, this->columns[col].combined.at(run_idx + 1));
        this->columns[col].combined.remove(run_idx);
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

        // Update pref_beg
        this->columns[col].pref_samples_beg.remove(run_idx);
        this->columns[col].pref_samples_beg.remove(run_idx);

        // Update pref_end
        this->columns[col].pref_samples_end.remove(run_idx);
        this->columns[col].pref_samples_end.remove(run_idx - 1);

        // Update div samples
        // Remove the current run and next run's div samples (Think about how deletion changes the div values)
        this->columns[col].div_samples_beg.remove(run_idx);
        this->columns[col].div_samples_beg.remove(run_idx);

        // Unset phi structures for hap_before and hap_after
        // For hap_before
        auto col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
        this->phi.phi_inv_supp[hap_before].remove(col_rank);
        this->phi.phi_inv_vec[hap_before].set(col, false);
        // For hap_after
        col_rank = this->phi.phi_vec[hap_after].rank1(col);
        this->phi.phi_supp[hap_after].remove(col_rank);
        this->phi.phi_supp_lcp[hap_after].remove(col_rank);
        this->phi.phi_vec[hap_after].set(col, false);

        // Unset Phi for hap_id
        col_rank = this->phi.phi_vec[hap_id].rank1(col);
        if (this->phi.phi_supp[hap_id].size() > 1) {
            this->phi.phi_supp[hap_id].remove(col_rank);
            this->phi.phi_supp_lcp[hap_id].remove(col_rank);
        }
        this->phi.phi_vec[hap_id].set(col, false);

        // Unset Phi_inv for hap_id
        col_rank = this->phi.phi_inv_vec[hap_id].rank1(col);
        if (this->phi.phi_inv_supp[hap_id].size() > 1)
            this->phi.phi_inv_supp[hap_id].remove(col_rank);
        this->phi.phi_inv_vec[hap_id].set(col, false);
    }

    void DeleteAtRunStart(const unsigned int col,
                          const unsigned int run_idx,
                          const unsigned int hap_id,
                          const bool allele) {
        unsigned int hap_before = 0;
        auto hap_after_opt = this->phi.phi_inv(hap_id, col);
        assert(hap_after_opt.has_value());
        unsigned int hap_after = hap_after_opt.value();
        unsigned int hap_after_div = 0;

        if (run_idx == 0) {
            hap_before = hap_after; // nothing above it
            hap_after_div = col;
        } else {
            hap_before = this->columns[col].pref_samples_end.at(run_idx - 1);
            hap_after_div = static_cast<int>(this->phi.plcp(hap_after, col));
        }

        // update run info spsi
        this->columns[col].combined.decrement(run_idx, 1);
        if (allele) {
            this->columns[col].ones.decrement(run_idx / 2, 1);
        } else {
            this->columns[col].zeros.decrement(run_idx / 2, 1);
            --this->columns[col].num_zeros;
        }

        // Update pref_beg with new head of run
        this->columns[col].pref_samples_beg.set(run_idx, hap_after);

        // Update phi for hap_id
        auto col_rank = this->phi.phi_vec[hap_id].rank1(col);
        if (this->phi.phi_supp[hap_id].size() > 1) {
            this->phi.phi_supp[hap_id].remove(col_rank);
            this->phi.phi_supp_lcp[hap_id].remove(col_rank);
        }
        this->phi.phi_vec[hap_id].set(col, false);

        // Update phi for new head of run
        col_rank = this->phi.phi_vec[hap_after].rank1(col);
        if (col != this->N) {
            this->phi.phi_supp[hap_after].insert(col_rank, hap_before);
        }
        this->phi.phi_vec[hap_after].set(col, true);

        // Update div for hap-after
        unsigned int curr_div = static_cast<int>(this->columns[col].div_samples_beg.at(run_idx));
        unsigned int new_div_val = max(curr_div, hap_after_div);
        this->columns[col].div_samples_beg.set(run_idx, new_div_val);
        if (col != this->N)
            this->phi.phi_supp_lcp[hap_after].insert(col_rank, new_div_val);

        // update phi for hap_before
        // checks if not at the top of a column
        if (hap_before != hap_after) {
            col_rank = this->phi.phi_inv_vec[hap_before].rank1(col);
            this->phi.phi_inv_supp[hap_before].set(col_rank, hap_after);
        }
    }

    void DeleteAtRunEnd(const unsigned int col,
                        const unsigned int run_idx,
                        const unsigned int hap_id,
                        const bool allele) {
        auto new_bottom_opt = this->phi.phi(hap_id, col);
        assert(new_bottom_opt.has_value());
        auto new_bottom = new_bottom_opt.value();
        unsigned int below_hap_id = new_bottom;
        if (run_idx + 1 < this->columns[col].pref_samples_beg.size()) {
            below_hap_id = this->columns[col].pref_samples_beg.at(run_idx + 1);
        }

        // Update pref_end
        this->columns[col].pref_samples_end.set(run_idx, new_bottom);

        // Update phi_inv for hap_id
        auto col_rank = this->phi.phi_inv_vec[hap_id].rank1(col);
        if (this->phi.phi_inv_supp[hap_id].size() > 1)
            this->phi.phi_inv_supp[hap_id].remove(col_rank);
        this->phi.phi_inv_vec[hap_id].set(col, false);

        // Update run_info
        this->columns[col].combined.decrement(run_idx, 1);
        if (allele) {
            this->columns[col].ones.decrement(run_idx / 2, 1);
        } else {
            this->columns[col].zeros.decrement(run_idx / 2, 1);
            --this->columns[col].num_zeros;
        }

        // Update phi_inv for new bottom and the haplotype below hap_id
        col_rank = this->phi.phi_inv_vec[new_bottom].rank1(col);
        if (col != this->N)
            this->phi.phi_inv_supp[new_bottom].insert(col_rank, below_hap_id);
        this->phi.phi_inv_vec[new_bottom].set(col, true);

        // Check if it's not at the bottom of the column
        if (below_hap_id != new_bottom) {
            col_rank = this->phi.phi_vec[below_hap_id].rank1(col);
            this->phi.phi_supp[below_hap_id].set(col_rank, new_bottom);

            // Update div value for the hap below
            unsigned int hap_id_div = this->phi.plcp(hap_id, col);
            unsigned int below_hap_id_div = this->columns[col].div_samples_beg.at(run_idx + 1);
            this->columns[col].div_samples_beg.set(run_idx + 1, max(hap_id_div, below_hap_id_div));
        }
    }

    void Delete(const unsigned int col, const unsigned int idx, const unsigned int hap_id, const bool allele) {
        unsigned int run_idx = get_run_idx(col, idx);
        if (isRunWithSingleElement(col, run_idx)) { // delete entire run
            DeleteSingleRun(col, idx, hap_id, allele);
        } else if (isRunStart(col, run_idx, hap_id)) { // deleting at head of a run
            DeleteAtRunStart(col, run_idx, hap_id, allele);
        } else if (isRunEnd(col, run_idx, hap_id)) { // deleting at end of a run
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

    void DeleteBottomHaplotype(const unsigned int hap_id) {
        // Handle deletion for the last haplotype
        vector<pair<unsigned int, unsigned int>> haplotype_info(this->N + 1); // {index, hapid}
        vector<bool> alleles;
        haplotype_info[0].first = hap_id;
        haplotype_info[0].second = hap_id;
        /*
        * find where each the target haplotype is mapped to in each column
        * Also find the allele values of the haplotype to be deleted
        */
        for (unsigned int col = 0; col < this->N; ++col) {
            unsigned int run_idx = get_run_idx(col, haplotype_info[col].first);
            const bool allele = this->get_run_val(col, run_idx);
            haplotype_info[col + 1] = w_mod(haplotype_info[col].first, col, allele, haplotype_info[col].second);
            alleles.push_back(allele);
        }
        assert(haplotype_info.size() == this->N + 1);
        assert(alleles.size() == this->N);

        // Perform deletion from each column
        for (unsigned int col = 0; col < this->N; ++col) {
            Delete(col, haplotype_info[col].first, haplotype_info[col].second, alleles[col]);
        }
        // Delete at the N-th column
        Delete(this->N, haplotype_info[this->N].first, haplotype_info[this->N].second, alleles[this->N - 1]);

        // checks
        for (unsigned int i = 0; i <= this->N; ++i) {
            assert(static_cast<bool>(this->phi.phi_vec[hap_id].at(i)) == false);
            assert(static_cast<bool>(this->phi.phi_inv_vec[hap_id].at(i)) == false);
        }
        assert(this->phi.phi_supp[hap_id].size() == 1);
        assert(this->phi.phi_inv_supp[hap_id].size() == 1);
        assert(this->phi.phi_supp_lcp[hap_id].size() == 1);

        // Handle for last column (Nth) update
        auto last_above = this->phi.phi_supp[hap_id].at(0);
        auto last_below = this->phi.phi_inv_supp[hap_id].at(0);
        // if removed haplotype not at the top in N-th col
        if (haplotype_info[this->N].first != 0) {
            assert(last_above != hap_id);
            // not at the bottom of the coln
            if (haplotype_info[this->N].first < (this->M - 1)) {
                assert(last_below != hap_id);
                this->phi.phi_inv_supp[last_above].set(this->phi.phi_inv_supp[last_above].size() - 1, last_below);
            } else {
                this->phi.phi_inv_supp[last_above].set(this->phi.phi_inv_supp[last_above].size() - 1, last_above);
            }
        }

        // if removed hap not at the bottom of the coln
        // update for the haplotype below it
        if (haplotype_info[this->N].first != (this->M - 1)) {
            assert(last_below != hap_id);
            if (haplotype_info[this->N].first > 0) {
                assert(last_above != hap_id);
                unsigned int max_div =
                    max(this->phi.phi_supp_lcp[last_below].at(this->phi.phi_supp_lcp[last_below].size() - 1),
                        this->phi.phi_supp_lcp[hap_id].at(this->phi.phi_supp_lcp[hap_id].size() - 1));
                this->phi.phi_supp[last_below].set(this->phi.phi_supp[last_below].size() - 1, last_above);
                this->phi.phi_supp_lcp[last_below].set(this->phi.phi_supp[last_below].size() - 1, max_div);
            } else {
                this->phi.phi_supp[last_below].set(this->phi.phi_supp[last_below].size() - 1, last_below);
                this->phi.phi_supp_lcp[last_below].set(this->phi.phi_supp[last_below].size() - 1, this->N);
            }
        }

        // Clean up phi supports
        this->phi.phi_vec.pop_back();
        this->phi.phi_supp.pop_back();
        this->phi.phi_inv_vec.pop_back();
        this->phi.phi_inv_supp.pop_back();
        this->phi.phi_supp_lcp.pop_back();
        --this->M;
        --this->phi.total_haplotypes;
    }

    void DeleteLastHaplotype(const unsigned int hap_id) {
        assert(hap_id == 0);
        for (unsigned int col = 0; col <= this->N; ++col) {
            // Update Run info
            this->columns[col].combined.remove(hap_id);
            assert(columns[col].combined.size() == 0);
            if (this->columns[col].zeros.size() == 0) {
                this->columns[col].ones.remove(hap_id);
                assert(this->columns[col].ones.size() == 0);
            } else {
                this->columns[col].zeros.remove(hap_id);
                assert(this->columns[col].zeros.size() == 0);
            }

            // Update Pref Beg/End info
            this->columns[col].pref_samples_beg.remove(0);
            assert(columns[col].pref_samples_beg.size() == 0);
            this->columns[col].pref_samples_end.remove(0);
            assert(columns[col].pref_samples_end.size() == 0);
            this->columns[col].div_samples_beg.remove(0);
            assert(columns[col].div_samples_beg.size() == 0);

            this->phi.phi_vec[hap_id].remove(0);
            this->phi.phi_inv_vec[hap_id].remove(0);
        }
        for (unsigned int i = 0; i <= this->N; ++i) {
            this->columns.pop_back();
        }
        assert(this->columns.size() == 0);
        assert(this->phi.phi_vec[hap_id].size() == 0);
        assert(this->phi.phi_inv_vec[hap_id].size() == 0);
        this->phi.phi_vec.pop_back();
        this->phi.phi_inv_vec.pop_back();

        // Update Phi Stuff
        while (this->phi.phi_supp[hap_id].size() > 0) {
            this->phi.phi_supp[hap_id].remove(0);
        }
        assert(this->phi.phi_supp[hap_id].size() == 0);
        while (this->phi.phi_inv_supp[hap_id].size() > 0) {
            this->phi.phi_inv_supp[hap_id].remove(0);
        }
        assert(this->phi.phi_inv_supp[hap_id].size() == 0);
        this->phi.phi_supp.pop_back();
        this->phi.phi_inv_supp.pop_back();
        --this->M;
        --this->phi.total_haplotypes;
    }

    /*
     * Input: haplotype index which will also be haplotype id for the first column (index is 0-based)
     * Result: Deletes that haplotype's allele from all columns
     */
    void DeleteSingleHaplotype(const unsigned int hap_id) {
        if (this->M == 0) {
            cerr << "Panel is empty. Nothing to delete! \n";
            return;
        }
        if (hap_id >= this->M) {
            cerr << "Invalid haplotype to delete! Hap ID should be in the range [0, " << this->M << ").\n";
            return;
        }
        if (this->M == 1) {
            DeleteLastHaplotype(hap_id);
            return;
        }

        // If deletion from the bottom
        if (hap_id == this->M - 1) {
            DeleteBottomHaplotype(hap_id);
        } else {
            // Get the last haplotype's alleles
            vector<pair<unsigned int, unsigned int>> last_hap_info(this->N + 1); // {index, hapid}
            vector<bool> last_hap_alleles;
            last_hap_info[0].first = this->M - 1;
            last_hap_info[0].second = this->M - 1;
            for (unsigned int col = 0; col < this->N; ++col) {
                unsigned int run_idx = get_run_idx(col, last_hap_info[col].first);
                const bool allele = this->get_run_val(col, run_idx);
                last_hap_info[col + 1] = w_mod(last_hap_info[col].first, col, allele, last_hap_info[col].second);
                last_hap_alleles.push_back(allele);
            }
            assert(last_hap_alleles.size() == this->N);
            last_hap_info.clear();
            /* Delete this bottom haplotype */
            DeleteBottomHaplotype(this->M - 1);

            /*
             * Delete The Target Haplotype
             * */
            vector<pair<unsigned int, unsigned int>> target_hap_info(this->N + 1); // {index, hapid}
            vector<bool> target_hap_alleles;
            target_hap_info[0].first = hap_id;
            target_hap_info[0].second = hap_id;
            for (unsigned int col = 0; col < this->N; ++col) {
                unsigned int run_idx = get_run_idx(col, target_hap_info[col].first);
                const bool allele = this->get_run_val(col, run_idx);
                target_hap_info[col + 1] = w_mod(target_hap_info[col].first, col, allele, target_hap_info[col].second);
                target_hap_alleles.push_back(allele);
            }
            assert(target_hap_alleles.size() == this->N);
            // Perform deletion from each column
            for (unsigned int col = 0; col < this->N; ++col) {
                Delete(col, target_hap_info[col].first, target_hap_info[col].second, target_hap_alleles[col]);
            }
            // Delete at the N-th column
            Delete(this->N,
                   target_hap_info[this->N].first,
                   target_hap_info[this->N].second,
                   target_hap_alleles[this->N - 1]);
            // checks
            for (unsigned int i = 0; i <= this->N; ++i) {
                assert(static_cast<bool>(this->phi.phi_vec[hap_id].at(i)) == false);
                assert(static_cast<bool>(this->phi.phi_inv_vec[hap_id].at(i)) == false);
            }
            assert(this->phi.phi_supp[hap_id].size() == 1);
            assert(this->phi.phi_inv_supp[hap_id].size() == 1);
            // Should this be true? apparently not, need to look further into
            // where supp_lcp is being updated when deletion is perfomred
            assert(this->phi.phi_supp_lcp[hap_id].size() == 1);

            // Handle for last column (Nth) update
            auto last_above = this->phi.phi_supp[hap_id].at(0);
            auto last_below = this->phi.phi_inv_supp[hap_id].at(0);

            // if removed haplotype not at the top in N-th col
            // Update phi_inv for hap above the removed haplotype
            if (target_hap_info[this->N].first != 0) {
                assert(last_above != hap_id);
                // not at the bottom of the coln
                if (target_hap_info[this->N].first < (this->M - 1)) {
                    assert(last_below != hap_id);
                    this->phi.phi_inv_supp[last_above].set(this->phi.phi_inv_supp[last_above].size() - 1, last_below);
                } else {
                    this->phi.phi_inv_supp[last_above].set(this->phi.phi_inv_supp[last_above].size() - 1, last_above);
                }
            }

            // if removed hap not at the bottom of the coln
            // Update phi for hap below the removed haplotype
            if (target_hap_info[this->N].first != (this->M - 1)) {
                assert(last_below != hap_id);
                if (target_hap_info[this->N].first > 0) {
                    assert(last_above != hap_id);
                    unsigned int max_div =
                        max(this->phi.phi_supp_lcp[last_below].at(this->phi.phi_supp_lcp[last_below].size() - 1),
                            this->phi.phi_supp_lcp[hap_id].at(this->phi.phi_supp_lcp[hap_id].size() - 1));
                    this->phi.phi_supp[last_below].set(this->phi.phi_supp[last_below].size() - 1, last_above);
                    this->phi.phi_supp_lcp[last_below].set(this->phi.phi_supp[last_below].size() - 1, max_div);
                } else {
                    this->phi.phi_supp[last_below].set(this->phi.phi_supp[last_below].size() - 1, last_below);
                    this->phi.phi_supp_lcp[last_below].set(this->phi.phi_supp[last_below].size() - 1, this->N);
                }
            }
            this->phi.phi_supp[hap_id].remove(0);
            this->phi.phi_inv_supp[hap_id].remove(0);
            // clear all elements from phi_supp_lcp
            this->phi.phi_supp_lcp[hap_id].remove(0);
            --this->M;
            --this->phi.total_haplotypes;

            /*
             * Insert the previous deleted last haplotype as hap_id
             */
            vector<pair<unsigned int, unsigned int>> insertion_indices(this->N + 1); // {index, hapid}
            insertion_indices[0].first = hap_id;
            if (this->M == 0) {
                insertion_indices[0].second = UINT_MAX;
            } else {
                insertion_indices[0].second = hap_id + 1;
            }
            for (unsigned int col = 0; col < this->N; ++col) {
                insertion_indices[col + 1] =
                    w_mod(insertion_indices[col].first, col, last_hap_alleles[col], insertion_indices[col].second);
            }

            vector<unsigned int> temp_div_seqn(this->N + 1, 0);
            vector<unsigned int> temp_div_below_seqn(this->N + 1, 0);
            CalculateDivergenceVal(insertion_indices, last_hap_alleles, temp_div_seqn, temp_div_below_seqn);
            packed_spsi temp_supp;
            packed_spsi temp_inv_supp;
            packed_spsi temp_div_supp;
            for (unsigned int col = 0; col <= this->N; ++col) {
                if (col == this->N) {
                    Insert(col,
                           insertion_indices[col].first,
                           insertion_indices[col].second,
                           hap_id,
                           last_hap_alleles[col - 1],
                           temp_supp,
                           temp_inv_supp,
                           temp_div_supp,
                           temp_div_seqn,
                           temp_div_below_seqn);

                } else {
                    Insert(col,
                           insertion_indices[col].first,
                           insertion_indices[col].second,
                           hap_id,
                           last_hap_alleles[col],
                           temp_supp,
                           temp_inv_supp,
                           temp_div_supp,
                           temp_div_seqn,
                           temp_div_below_seqn);
                }
            }

            // Final Updates
            this->phi.phi_supp[hap_id] = temp_supp;
            this->phi.phi_inv_supp[hap_id] = temp_inv_supp;
            this->phi.phi_supp_lcp[hap_id] = temp_div_supp;
            ++this->M;
            ++this->phi.total_haplotypes;
        }
    }

    /*
     * Long match query algorithm
     */
    std::vector<std::tuple<int, int, int>> compute_long_match(vector<bool> &query,
                                                              const unsigned int length,
                                                              bool verbose = false) {
        // compute the match iff |query| is equal to the width of the panel
        if (query.size() != this->N) {
            std::cout << query.size() << " != " << this->N
                      << "\n";
            exit(1);
        }

        // virtual insertion
        std::vector<std::pair<unsigned int, unsigned int>> virtual_insert_indices(this->N + 1);
        std::vector<unsigned int> zstart(this->N + 1);
        std::vector<unsigned int> bstart(this->N + 1);

        // update sequence below query in each column in
        virtual_insert_indices[0].first = this->M; // stores haplotype index
        virtual_insert_indices[0].second = UINT_MAX; // stores haplotype ID
        for (unsigned int k = 0; k < this->N; ++k) {
            virtual_insert_indices[k + 1] =
                w_mod(virtual_insert_indices[k].first, k, query[k], virtual_insert_indices[k].second);
        }

        // update divergence values
        CalculateDivergenceVal(virtual_insert_indices, query, zstart, bstart);

        // query algorithm
        std::pair<unsigned int, unsigned int> f_curr_pair = virtual_insert_indices[0];
        std::pair<unsigned int, unsigned int> g_curr_pair = virtual_insert_indices[0];
        std::vector<unsigned int> dZ(this->M, 0); // store div values

        // stores the matches found
        std::vector<std::tuple<int, int, int>> matches;
        for (unsigned int k = 0; k < this->N; ++k) {
            bool query_opposite = !(query[k]);
            auto f_opp_pair = w_mod(f_curr_pair.first, k, query_opposite, f_curr_pair.second);
            auto g_opp_pair = w_mod(g_curr_pair.first, k, query_opposite, g_curr_pair.second);

            auto f_prime_pair = w_mod(f_curr_pair.first, k, query[k], f_curr_pair.second);
            auto g_prime_pair = w_mod(g_curr_pair.first, k, query[k], g_curr_pair.second);

            while (f_opp_pair.first != g_opp_pair.first) {
                // store matches
                matches.emplace_back(f_opp_pair.second, dZ[f_opp_pair.second], k);

                // equivalent to : f_temp  = f_temp.below
                auto new_f_temp_hapid = this->phi.phi_inv(f_opp_pair.second, k + 1);
                f_opp_pair.first++;
                if (new_f_temp_hapid.has_value()) {
                    f_opp_pair.second = new_f_temp_hapid.value();
                } else {
                    assert(f_opp_pair.first == g_opp_pair.first);
//                    if (verbose)
//                        std::cout << "Report matches: Bottom boundary reached \n";
                    break;
                }
            }

            // if empty block
            assert((f_prime_pair.first == g_prime_pair.first) == (f_prime_pair.second == g_prime_pair.second));
            if (f_prime_pair.first == g_prime_pair.first) {
                if (k + 1 - zstart[k + 1] == length) {
                    assert(f_prime_pair.first > 0);
                    auto above_val = this->phi.phi(f_prime_pair.second, k + 1);
                    if (above_val.has_value()) {
                        f_prime_pair.second = above_val.value();
                        assert(f_prime_pair.first > 0);
                        f_prime_pair.first = std::max((int) f_prime_pair.first - 1, 0);
                    } else {
                        std::cout << "Error in finding phi for hapID : " << f_prime_pair.second << " in column " << k
                                  << "\n";
                        exit(EXIT_FAILURE);
                    }
                    dZ[f_prime_pair.second] = k + 1 - length;
                }

                if (k + 1 - bstart[k + 1] == length) {
                    dZ[g_prime_pair.second] = k + 1 - length;
                    auto below_val = this->phi.phi_inv(g_prime_pair.second, k + 1);
                    if (below_val.has_value()) {
                        g_prime_pair.second = below_val.value();
                        assert(g_prime_pair.first < this->M);
                        g_prime_pair.first = std::max((int) g_prime_pair.first + 1, 0);
                    } else {
                        std::cout << "Error in finding phi for hapID : " << g_prime_pair.second << " in column " << k
                                  << "\n";
                        exit(EXIT_FAILURE);
                    }
                }
            }

            // expand boundaries of the block
            if (f_prime_pair.first != g_prime_pair.first) {
                while (this->phi.plcp(f_prime_pair.second, k + 1) <= k + 1 - length) {
                    auto above_val = this->phi.phi(f_prime_pair.second, k + 1);
                    if (above_val.has_value()) {
                        f_prime_pair.second = above_val.value();
                        f_prime_pair.first = std::max((int) f_prime_pair.first - 1, 0);
                    } else {
                        std::cout << "Error in finding phi for hapID : " << f_prime_pair.second << " in column " << k
                                  << "\n";
                        exit(EXIT_FAILURE);
                        if (verbose)
                            std::cout << "Expand boundary: Top boundary reached \n";
                        break;
                    }
                    dZ[f_prime_pair.second] = k + 1 - length;
                }
                while (this->phi.plcp(g_prime_pair.second, k + 1) <= k + 1 - length) {
                    dZ[g_prime_pair.second] = k + 1 - length;
                    auto below_val = this->phi.phi_inv(g_prime_pair.second, k + 1);
                    if (below_val.has_value()) {
                        g_prime_pair.second = below_val.value();
                        assert(g_prime_pair.first < this->M);
                        g_prime_pair.first = std::min((int) g_prime_pair.first + 1, (int) this->M);
                    } else {
                        if (verbose)
                            std::cout << "Expand boundary: Bottom Boundary reached: \n";
                        g_prime_pair.second = -1;
                        g_prime_pair.first = (int) this->M;
                        break;
                    }
                }
            }
            f_curr_pair = f_prime_pair;
            g_curr_pair = g_prime_pair;
        }
        while (f_curr_pair.first != g_curr_pair.first) {
            // store matches
            matches.emplace_back(f_curr_pair.second, dZ[f_curr_pair.second], this->N);
            // equivalent to : f_temp  = f_temp.below
            auto new_f_temp_hapid = this->phi.phi_inv(f_curr_pair.second, this->N);
            f_curr_pair.first++;
            if (new_f_temp_hapid.has_value()) {
                f_curr_pair.second = new_f_temp_hapid.value();
            } else {
                assert(f_curr_pair.first == this->M && f_curr_pair.first == g_curr_pair.first);
                if (verbose)
                    std::cout << "Report matches: Bottom boundary reached \n";
                break;
            }
        }
        return matches;
    }

    // TODO: Test this!!!
    std::vector<std::tuple<int, int, int>> compute_set_maximal_match(vector<bool> &query) {
        // compute the match iff |query| is equal to the width of the panel
        if (query.size() != this->N) {
            std::cout << query.size() << " != " << this->N
                      << "\n";
            exit(1);
        }

        std::vector<std::pair<unsigned int, unsigned int>> virtual_insert_indices(this->N + 1);
        std::vector<unsigned int> zstart(this->N + 1);
        std::vector<unsigned int> bstart(this->N + 1);

        // update sequence below query in each column in
        virtual_insert_indices[0].first = this->M; // stores haplotype index
        virtual_insert_indices[0].second = UINT_MAX; // stores haplotype ID
        for (unsigned int k = 0; k < this->N; ++k) {
            virtual_insert_indices[k + 1] =
                w_mod(virtual_insert_indices[k].first, k, query[k], virtual_insert_indices[k].second);
        }

        // update divergence values
        CalculateDivergenceVal(virtual_insert_indices, query, zstart, bstart);
        auto outputMatches = [this](unsigned site,
                                    unsigned hapId,
                                    unsigned matchStart,
                                    std::vector<std::tuple<int, int, int>> &matches) {
            //output matches to hapId that match on [matchStart, site)
            matches.emplace_back(hapId, matchStart, site);
            unsigned topOfBlock = hapId;
            //output matches above
            while (this->phi.plcp(topOfBlock, site) <= matchStart) {
                auto above_val = this->phi.phi(topOfBlock, site);
                assert(above_val.has_value());
                if (above_val.value() == topOfBlock)
                    //std::cout << "above_val == topOfBlock == " << above_val.value() << std::endl;
                    topOfBlock = above_val.value();
                matches.emplace_back(topOfBlock, matchStart, site);
            }
            //std::cout << "expanding block below" << std::endl;
            //output matches below
            unsigned botOfBlock = hapId;
            auto belowBlock = this->phi.phi_inv(hapId, site);
            //std::cout << "starting while loop" << std::endl;
            while (botOfBlock != belowBlock.value() && this->phi.plcp(belowBlock.value(), site) <= matchStart) {
                matches.emplace_back(belowBlock.value(), matchStart, site);
                botOfBlock = belowBlock.value();
                belowBlock = this->phi.phi_inv(botOfBlock, site);
            }
            //std::cout << "ending outputMatches" << std::endl;
        };

        unsigned prevDiv = 0, newDiv;
        std::vector<std::tuple<int, int, int>> matches;
        assert(prevDiv == std::min(zstart[0], bstart[0]));
        for (unsigned k = 1; k <= this->N; ++k, prevDiv = newDiv) {
            newDiv = std::min(zstart[k], bstart[k]);
            if (newDiv == prevDiv || prevDiv == k - 1)
                continue;
            //std::cout << "computing hapId" << std::endl;
            unsigned hapId = virtual_insert_indices[k - 1].second;
            if (zstart[k - 1] < bstart[k - 1]) {
                hapId = (hapId == UINT_MAX) ? this->columns[k - 1].pref_samples_end[
                    this->columns[k - 1].pref_samples_end.size() - 1] : this->phi.phi(hapId, k - 1).value();
            }
            //std::cout << "finished computing hapId" << std::endl;
            outputMatches(k - 1, hapId, prevDiv, matches);
        }
        //std::cout << "outside of while loop" << std::endl;
        if (prevDiv != this->N) {
            unsigned k = this->N + 1;
            unsigned hapId = (zstart[k - 1] < bstart[k - 1]) ?
                             this->phi.phi(virtual_insert_indices[k - 1].second, k - 1).value() :
                             virtual_insert_indices[k - 1].second;
            outputMatches(k - 1, hapId, prevDiv, matches);
        }
        //std::cout << "returning from smem query" << std::endl;
        return matches;
    }

    // set max match query
    void smem_query(string &query_vcf, string &out) {
        std::ofstream out_match(out);
        vector<vector<bool>> queries;
        ReadQueryFile(query_vcf.c_str(), queries);
        unsigned int n_queries = queries.size();
        std::vector<std::vector<std::tuple<int, int, int>>> matches_vec(n_queries);
        std::cout << "starting smem queries" << std::endl;
        clock_t START = clock();
        for (unsigned int i = 0; i < n_queries; i++) {
            clock_t st = clock();
            matches_vec[i] = this->compute_set_maximal_match(queries[i]);
            auto end = (float) (clock() - st) / CLOCKS_PER_SEC;
            out_match << end << "\n";
            cout << "Completed query " << i << " in " << end << " seconds" << std::endl;
        }
        auto query_time = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Time to query = " << query_time << " s.\n";
        std::cout << "ending smem queries" << std::endl;

        // report long matches
        START = clock();
        std::cout << "Outputting matches...\n";
        out_match << "MATCH\tQueryHapID\tRefHapID\tStart\tEnd(exclusive)\tLength\n";
        for (unsigned int i = 0; i < queries.size(); i++) {
            for (unsigned int j = 0;
                 j < matches_vec[i].size(); j++) {
                auto ref_hap_ID =
                    std::get<0>(
                        matches_vec[i][j]);
                auto start_pos =
                    std::get<1>(
                        matches_vec[i][j]);
                auto end =
                    std::get<2>(
                        matches_vec[i][j]);
                out_match << "MATCH\t" << i << "\t" << ref_hap_ID << "\t"
                          << start_pos << "\t" << end << "\t"
                          << end - start_pos << "\n";
            }
        }
        auto time_write_matches = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Time to output matches = " << time_write_matches << " s.\n";
        out_match.close();
    }

    void long_match_query(string &query_vcf, string &out, unsigned int length = 0,
                          bool verbose = false) {
        std::ofstream out_match(out);

        vector<vector<bool>> queries;
        ReadQueryFile(query_vcf.c_str(), queries);
        unsigned int n_queries = queries.size();
        std::vector<std::vector<std::tuple<int, int, int>>> matches_vec(n_queries);
        clock_t START = clock();
#pragma omp parallel for default(none) \
                       shared(queries, matches_vec, n_queries, length, verbose, std::cout)
        for (unsigned int i = 0; i < n_queries; i++) {
            //clock_t st = clock();
            matches_vec[i] = this->compute_long_match(queries[i], length, verbose);
            //auto end = (float)(clock() - st)/CLOCKS_PER_SEC;
            //out_match << end << "\n";
        }
        auto query_time = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Time to query = " << query_time << " s.\n";

        // report long matches
        START = clock();
        std::cout << "Outputting matches...\n";
        out_match << "MATCH\tQueryHapID\tRefHapID\tStart\tEnd(exclusive)\tLength\n";
        for (unsigned int i = 0; i < queries.size(); i++) {
            for (unsigned int j = 0;
                 j < matches_vec[i].size(); j++) {
                auto ref_hap_ID =
                    std::get<0>(
                        matches_vec[i][j]);
                auto start_pos =
                    std::get<1>(
                        matches_vec[i][j]);
                auto end =
                    std::get<2>(
                        matches_vec[i][j]);
                out_match << "MATCH\t" << i << "\t" << ref_hap_ID << "\t"
                          << start_pos << "\t" << end << "\t"
                          << end - start_pos << "\n";
            }
        }
        auto time_write_matches = (float) (clock() - START) / CLOCKS_PER_SEC;
        std::cout << "Time to output matches = " << time_write_matches << " s.\n";
        out_match.close();
    }

    unsigned long long get_memory_usage_bytes() {
        unsigned long long pa_da_spsi_size_bits = 0, run_info_spsi_size_bits = 0;
        unsigned long long others = 0;
        for (unsigned int col = 0; col <= this->N; ++col) {
            run_info_spsi_size_bits += this->columns[col].combined.bit_size();
            run_info_spsi_size_bits += this->columns[col].zeros.bit_size();
            run_info_spsi_size_bits += this->columns[col].ones.bit_size();
            pa_da_spsi_size_bits += this->columns[col].pref_samples_beg.bit_size();
            pa_da_spsi_size_bits += this->columns[col].pref_samples_end.bit_size();
            pa_da_spsi_size_bits += this->columns[col].div_samples_beg.bit_size();
        }

        others += sizeof(bool) * (this->N + 1); // start_with_zero
        others += sizeof(unsigned int) * (this->N + 1); // num_zeros
//    others /= 8;

        unsigned long long pa_da_spsi_size_bytes = pa_da_spsi_size_bits / 8;
        unsigned long long run_info_spsi_size_bytes = run_info_spsi_size_bits / 8;
        unsigned long long all_column_spsi_memory = pa_da_spsi_size_bytes + run_info_spsi_size_bytes + others;
        unsigned long long phi_memory = this->phi.size_in_bytes(false);
        return all_column_spsi_memory + phi_memory;
    }

    void PrintMemoryUsage(bool verbose = false) {
        unsigned long long pa_da_spsi_size_bits = 0, run_info_spsi_size_bits = 0;
        unsigned long long others = 0;
        for (unsigned int col = 0; col <= this->N; ++col) {
//      assert(this->columns[col].combined.size());
//       assert(this->columns[col].zeros.size());
//       assert(this->columns[col].ones.size());
//      assert(this->columns[col].pref_samples_beg.size());
//      assert(this->columns[col].pref_samples_end.size());
//      assert(this->columns[col].div_samples_beg.size());
            run_info_spsi_size_bits += this->columns[col].combined.bit_size();
            run_info_spsi_size_bits += this->columns[col].zeros.bit_size();
            run_info_spsi_size_bits += this->columns[col].ones.bit_size();
            pa_da_spsi_size_bits += this->columns[col].pref_samples_beg.bit_size();
            pa_da_spsi_size_bits += this->columns[col].pref_samples_end.bit_size();
            pa_da_spsi_size_bits += this->columns[col].div_samples_beg.bit_size();
        }

        others += sizeof(bool) * (this->N + 1); // start_with_zero
        others += sizeof(unsigned int) * (this->N + 1); // num_zeros
//    others /= 8;

        unsigned long long pa_da_spsi_size_bytes = pa_da_spsi_size_bits / 8;
        unsigned long long run_info_spsi_size_bytes = run_info_spsi_size_bits / 8;
        unsigned long long all_column_spsi_memory = pa_da_spsi_size_bytes + run_info_spsi_size_bytes + others;
        unsigned long long phi_memory = this->phi.size_in_bytes(verbose);
        std::cout << "\n";
        std::cout << "##### Memory consumption #####" << std::endl;
        std::cout << "Memory of dcPBWT pa/da spsis : " << pa_da_spsi_size_bytes << " bytes." << std::endl;
        std::cout << "Memory of dcPBWT run_info (combined + ones + zeros) spsis : " << run_info_spsi_size_bytes
                  << " bytes." << std::endl;
        std::cout << "Memory of dcPBWT others (start_with_zero+num_zeros): " << others << " bytes." << std::endl;
        std::cout << "Total Memory of dcPBWT columns spsis : " << all_column_spsi_memory << " bytes." << std::endl;
        std::cout << "Memory of dcPBWT phi : " << phi_memory << " bytes." << std::endl;
        std::cout << "Total memory of dcPBWT: " << phi_memory + all_column_spsi_memory << " bytes." << std::endl;
        std::cout << std::endl;
    }

    double get_avg_runs() {
        double total_runs = 0;
        for (unsigned int col = 0; col < this->N; ++col)
            total_runs += this->columns[col].combined.size();
        double avg_runs = total_runs / this->N;
        return avg_runs;
    }

    void load(std::istream &in) {
        in.read((char *) &this->N, sizeof(this->N));
        in.read((char *) &this->M, sizeof(this->M));
        this->columns = std::vector<dcpbwt_column>(this->N + 1);
        for (unsigned int col = 0; col < this->columns.size(); ++col) {
            this->columns[col].combined.load(in);
            this->columns[col].ones.load(in);
            this->columns[col].zeros.load(in);
            this->columns[col].pref_samples_beg.load(in);
            this->columns[col].pref_samples_end.load(in);
            this->columns[col].div_samples_beg.load(in);
            in.read((char *) &this->columns[col].start_with_zero, sizeof(bool));
            in.read((char *) &this->columns[col].num_zeros, sizeof(unsigned int));
        }
        this->phi.load(in);
    }

    size_t serialize(std::ostream &out) const {
        size_t written_bytes = 0;
        out.write((char *) &this->N, sizeof(this->N));
        written_bytes += sizeof(this->N);
        out.write((char *) &this->M, sizeof(this->M));
        written_bytes += sizeof(this->M);

        for (unsigned int col = 0; col < this->columns.size(); ++col) {
            written_bytes += this->columns[col].combined.serialize(out);
            written_bytes += this->columns[col].ones.serialize(out);
            written_bytes += this->columns[col].zeros.serialize(out);
            written_bytes += this->columns[col].pref_samples_beg.serialize(out);
            written_bytes += this->columns[col].pref_samples_end.serialize(out);
            written_bytes += this->columns[col].div_samples_beg.serialize(out);
            out.write((char *) &this->columns[col].start_with_zero, sizeof(this->columns[col].start_with_zero));
            written_bytes += sizeof(bool);
            out.write((char *) &this->columns[col].num_zeros, sizeof(this->columns[col].num_zeros));
            written_bytes += sizeof(unsigned int);
        }
        written_bytes += phi.serialize(out);
        return written_bytes;
    }

    // Build the reference panel using htslib
    void Build(const char *filename) {
        htsFile *fp = hts_open(filename, "rb");
        std::cout << "Building panel from: " << filename << std::endl;
        if (fp == nullptr) {
            std::cerr << "VCF file: " << filename << " not found!" << std::endl;
            return;
        }
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        bcf1_t *rec = bcf_init();
        this->M = bcf_hdr_nsamples(hdr) * 2;
        fp = hts_open(filename, "rb");
        hdr = bcf_hdr_read(fp);
        rec = bcf_init();

        // go through all sites
        vector<unsigned int> u, v;
        vector<unsigned int> freq;
        unsigned int total_runs = 0;
        unsigned int col = 0;
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

        // iterate each vcf record
        vector<bool> single_col;
        while (bcf_read(fp, hdr, rec) >= 0) {
            this->N++;
            single_col.clear();
            bcf_unpack(rec, BCF_UN_ALL);
            // read SAMPLE
            int32_t *gt_arr = nullptr, ngt_arr = 0;
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
                    single_col.push_back(allele_index == 1);
                }
            }
            // Build column
            packed_spsi temp_zeros;
            packed_spsi temp_ones;
            packed_spsi temp_combined;
            packed_spsi temp_sample_beg;
            packed_spsi temp_sample_end;
            packed_spsi temp_div;
            bool start_with_zero = false;
            assert(single_col.size() == this->M);
            p = col + 1;
            q = col + 1;
            for (unsigned int i = 0; i < this->M; ++i) {
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
            for (unsigned int i = 0; i < freq.size(); ++i) {
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
            for (unsigned int i = 0; i < temp_sample_beg.size(); ++i) {
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

            free(gt_arr);
        }

        // Handle for last column
        packed_spsi temp_zeros;
        packed_spsi temp_ones;
        packed_spsi temp_combined;
        packed_spsi temp_sample_beg;
        packed_spsi temp_sample_end;
        packed_spsi temp_div;
        bool start_with_zero = false;
        assert(single_col.size() == this->M);
        p = col + 1;
        q = col + 1;
        for (unsigned int i = 0; i < this->M; ++i) {
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
        temp_sample_end.push_back(prefix_arr[this->M - 1]);
        // edge case
        if (cnt > 0)
            freq.push_back(cnt);

        // populate dynamic data structures
        for (unsigned int i = 0; i < freq.size(); ++i) {
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
        for (unsigned int i = 0; i < temp_sample_beg.size(); ++i) {
            sites_where_sample_beg[temp_sample_beg.at(i)].push_back(col);
            sites_where_sample_end[temp_sample_end.at(i)].push_back(col);
        }
        // build each column
        dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                           (std::move(temp_sample_beg)), (std::move(temp_sample_end)), (std::move(temp_div)),
                           start_with_zero, u.size());
        columns.emplace_back(coln);
        assert(this->columns.size() == this->N + 1);
        this->phi = phi_ds(columns, M, N, sites_where_sample_beg, sites_where_sample_end, prefix_arr, div);
        assert(col == N);
        total_runs += freq.size();
        cout << "Avg runs in panel = " << static_cast<float>(total_runs) / this->N << "\n";
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        bcf_destroy(rec);
    }

    // Build the reference panel
    void BuildFromVCF(std::string &filename) {
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
                this->M += 2;
                // individual_ids.push_back(token);
            }
            // go through all sites
            vector<unsigned int> u{}, v{};
            vector<unsigned int> freq{};
            unsigned int total_runs = 0;
            unsigned int col = 0;
            int cnt = 1;
            bool prev_allele = false;
            // TODO: Possible optimization using sdsl int vec
            vector<unsigned int> prefix_arr(this->M, 0);
            std::iota(prefix_arr.begin(), prefix_arr.end(), 0);
            // TODO: Possible optimization using sdsl int vec
            vector<vector<unsigned int>> sites_where_sample_beg(this->M);
            vector<vector<unsigned int>> sites_where_sample_end(this->M);
            vector<unsigned int> div(this->M, 0);
            vector<unsigned int> div_u{}, div_v{};
            unsigned int p = col + 1, q = col + 1;

            // read through the line and store in single_col
            std::vector<bool> single_col;
            while (getline(inFile, line)) {
                std::istringstream iss(line);
                token = "";
                single_col.clear();
                for (unsigned int i = 0; i < (this->M / 2) + 9; ++i) {
                    iss >> token;
                    if (i < 9) {
                        continue;
                    }
                    single_col.push_back(static_cast<bool>(token[0] - '0'));
                    single_col.push_back(static_cast<bool>(token[2] - '0'));
                }

                packed_spsi temp_zeros{};
                packed_spsi temp_ones{};
                packed_spsi temp_combined{};
                packed_spsi temp_sample_beg{};
                packed_spsi temp_sample_end{};
                packed_spsi temp_div{};
                bool start_with_zero = false;
                assert(single_col.size() == this->M);
                p = col + 1;
                q = col + 1;
                for (unsigned int i = 0; i < this->M; ++i) {
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
                for (unsigned int i = 0; i < freq.size(); ++i) {
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

                for (unsigned int i = 0; i < temp_sample_beg.size(); ++i) {
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

            this->N = col;
            // Handle for last column boundary
            packed_spsi temp_zeros{};
            packed_spsi temp_ones{};
            packed_spsi temp_combined{};
            packed_spsi temp_sample_beg{};
            packed_spsi temp_sample_end{};
            packed_spsi temp_div{};
            bool start_with_zero = false;
            assert(single_col.size() == this->M);
            p = col + 1;
            q = col + 1;
            for (unsigned int i = 0; i < this->M; ++i) {
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
            temp_sample_end.push_back(prefix_arr[this->M - 1]);
            // edge case
            if (cnt > 0)
                freq.push_back(cnt);

            // populate dynamic data structures
            for (unsigned int i = 0; i < freq.size(); ++i) {
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
            for (unsigned int i = 0; i < temp_sample_beg.size(); ++i) {
                sites_where_sample_beg[temp_sample_beg.at(i)].push_back(col);
                sites_where_sample_end[temp_sample_end.at(i)].push_back(col);
            }
            // build each column
            dcpbwt_column coln((std::move(temp_zeros)), (std::move(temp_ones)), (std::move(temp_combined)),
                               (std::move(temp_sample_beg)), (std::move(temp_sample_end)), (std::move(temp_div)),
                               start_with_zero, u.size());
            columns.emplace_back(coln);
            assert(columns.size() == this->N + 1);

            // build phi data-structure
            this->phi =
                phi_ds(columns, this->M, this->N, sites_where_sample_beg, sites_where_sample_end, prefix_arr, div);
            assert(col == N);
            total_runs += freq.size();
            cout << "Avg runs = " << static_cast<float>(total_runs) / this->N << "\n";

            inFile.close();
        } else {
            std::cerr << "Couldn't find : " << filename << "\n";
            exit(1);
        }
    }
};
