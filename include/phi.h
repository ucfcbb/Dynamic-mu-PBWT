//
// Created by shakyap on 8/9/24.
//

#ifndef PHI_H
#define PHI_H

#include <iostream>
#include <optional>

#include "dcpbwt_column.h"
#include "dynamic/dynamic.hpp"

using namespace dyn;

class phi_ds{
private:
    unsigned int total_haplotypes{};
    unsigned int total_sites{};
public:
    /**
     * sparse dynamic sparse bitvectors for phi function
     */
    std::vector<suc_bv> phi_vec;

    /**
     * panel of dynamic sparse bitvectors for phi_inv function
     */
    std::vector<suc_bv> phi_inv_vec;


    /**
     * @brief compressed int vector for prefix samples used by phi function
     */
    std::vector<packed_spsi> phi_supp;

    /**
     * @brief compressed int vector for prefix samples used by phi_inv function
     */
    std::vector<packed_spsi> phi_inv_supp;

    /**
     * @brief default constructor
     */
    phi_ds() = default;

    /**
     * @brief default destructor
     */
    virtual ~phi_ds() = default;

    /**
     * @brief constructor of the phi/phi_inv support data structure
     * @param cols vector of the columns of the RLPBWT
     * @param panelbv random access data structure for the panel of RLPBWT
     * @param last_pref last prefix array of the PBWT
     * @param verbose bool for extra prints
     */
    phi_ds(vector<dcpbwt_column>& columns, unsigned int M,
                    unsigned int N,
                    // sdsl::int_vector<> &last_div,
                    std::vector<std::vector<unsigned int>> &site_where_samples_beg,
                    std::vector<std::vector<unsigned int>> &site_where_samples_end,
                    vector<int> &last_pref,
                    bool verbose = false) {
        // default value is the panel height
        this->total_haplotypes = M;
        this->total_sites = N;

        // initialize panels and vectors of the data structure
        this->phi_vec = std::vector<suc_bv>(total_haplotypes);
        this->phi_inv_vec = std::vector<suc_bv>(total_haplotypes);

        this->phi_supp = std::vector<packed_spsi>(total_haplotypes);
        this->phi_inv_supp = std::vector<packed_spsi>(total_haplotypes);

        // temporary vector for supports
        // std::vector<std::vector<unsigned int>> phi_supp_tmp(total_haplotypes);
        // std::vector<std::vector<unsigned int>> phi_inv_supp_tmp(total_haplotypes);

        // build sparse bitvector
        for (unsigned int hap = 0; hap < total_haplotypes; hap++) {
            suc_bv tmp_begin, tmp_end;
            // initialize both suc_bv's with 0s (i.e. false)
            for(int i = 0; i < total_sites + 1; ++i) {
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
        for (unsigned int col = 0; col < this->total_sites; col++) {
            for (unsigned int j = 0; j < columns[col].pref_samples_beg.size(); j++) {
                // use sample beg to compute phi panel
                // use sample_end to compute
                // support phi panel (if we are in the first run we use default
                // value)
                if (j == 0) {
                    this->phi_supp[columns[col].pref_samples_beg[j]].push_back(total_haplotypes);
                    // phi_supp_tmp_l[cols[i].sample_beg[j]].push_back(0);
                } else {
                    this->phi_supp[columns[col].pref_samples_beg[j]].push_back(columns[col].pref_samples_end[j - 1]);
                    // phi_supp_tmp[cols[i].sample_beg[j]].push_back(
                    //         cols[i].sample_end[j - 1]);
                    // phi_supp_tmp_l[cols[i].sample_beg[j]].push_back(
                    //         cols[i].sample_beg_lcp[j]);
                }
                // use sample end to compute phi_inv panel
                // use sample_beg to compute
                // support phi panel (if we are in the last run we use default
                // value)
                if (j == columns[col].pref_samples_beg.size() - 1) {
                    this->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
                            total_haplotypes);
                } else {
                    this->phi_inv_supp[columns[col].pref_samples_end[j]].push_back(
                            columns[col].pref_samples_beg[j+1]);
                }
            }
        }
        // use the last prefix array to compute the remain values for the
        // phi support data structure (with the same "rules" of the previous
        // case)
        for (unsigned int j = 0; j < this->phi_supp.size(); j++) {
            if (j == 0) {
                if ((this->phi_supp[last_pref[j]].size() == 0) ||
                    this->phi_supp[last_pref[j]][this->phi_supp[last_pref[j]].size() - 1] != this->total_haplotypes) {
                        this->phi_supp[last_pref[j]].push_back(total_haplotypes);
                    // phi_supp_tmp_l[last_pref[j]].push_back(0);
                }
            } else {
                if ((this->phi_supp[last_pref[j]].size() == 0) ||
                    this->phi_supp[last_pref[j]][this->phi_supp[last_pref[j]].size() - 1] != last_pref[j - 1]) {
                    this->phi_supp[last_pref[j]].push_back(last_pref[j - 1]);
                    // phi_supp_tmp_l[last_pref[j]].push_back(last_div[j]);
                }
            }
            if (j == this->phi_supp.size() - 1) {
                if (this->phi_inv_supp[last_pref[j]].size() == 0 ||
                    this->phi_inv_supp[last_pref[j]][this->phi_inv_supp[last_pref[j]].size() - 1] != total_haplotypes) {
                    this->phi_inv_supp[last_pref[j]].push_back(total_haplotypes);
                }
            } else {
                if (this->phi_inv_supp[last_pref[j]].size() == 0 ||
                    this->phi_inv_supp[last_pref[j]][this->phi_inv_supp[last_pref[j]].size() - 1] != last_pref[j + 1]) {
                    this->phi_inv_supp[last_pref[j]].push_back(last_pref[j + 1]);
                }
            }
        }
    }

    /**
     * @brief phi function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return previous prefix array value at current column (if exists)
     */
    std::optional<unsigned int> phi(unsigned int pref, unsigned int col) {
        // assert(col < this->total_sites);
        auto tmp_col = this->phi_vec[pref].rank1(col);
        if(tmp_col == this->phi_supp[pref].size()){
            tmp_col--;
        }
        auto res = static_cast<unsigned int>(this->phi_supp[pref].at(tmp_col));
        // auto res = this->phi_supp[pref].at(tmp_col);
        cout << res << "\n";
        if (res == this->total_haplotypes) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    /**
     * @brief phi_inv function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return next prefix array value at current column (if exists)
    */
    std::optional<unsigned int> phi_inv(unsigned int pref, unsigned int col) {
        auto tmp_col = this->phi_inv_vec[pref].rank1(col);
        if(tmp_col == this->phi_inv_supp[pref].size()){
            tmp_col--;
        }
        auto res = static_cast<unsigned int>(this->phi_inv_supp[pref][tmp_col]);
        if (res == this->total_haplotypes) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    // unsigned int plcp(unsigned int pref, unsigned int col) {
    //     if (col == 0 || !phi(pref, col).has_value()) {
    //         return 0;
    //     }
    //     auto tmp_col = this->phi_rank[pref](col);
    //     if(tmp_col == this->phi_supp[pref].size()){
    //         tmp_col--;
    //     }
    //     auto end_col = this->phi_select[pref](tmp_col + 1);
    //     auto tmp = static_cast<int>(this->phi_supp_l[pref][tmp_col]);
    //     auto plcp = tmp - (end_col - col);
    //     return plcp;
    // }

    /**
     * @brief function to obtain size in bytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    // unsigned long long size_in_bytes(bool verbose = false) {
    //     unsigned long long size = 0;
    //     for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
    //         size += sdsl::size_in_bytes(phi_vec[i]);
    //         size += sdsl::size_in_bytes(phi_inv_vec[i]);
    //         size += sdsl::size_in_bytes(phi_rank[i]);
    //         size += sdsl::size_in_bytes(phi_select[i]);
    //         size += sdsl::size_in_bytes(phi_inv_rank[i]);
    //         size += sdsl::size_in_bytes(phi_inv_select[i]);
    //         size += sdsl::size_in_bytes(phi_supp[i]);
    //         size += sdsl::size_in_bytes(phi_inv_supp[i]);
    //         size += sdsl::size_in_bytes(phi_supp_l[i]);
    //     }
    //     if (verbose) {
    //         std::cout << "phi support: " << size << " bytes\n";
    //     }
    //     return size;
    // }

    /**
     * @brief function to obtain size in megabytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    // double size_in_mega_bytes(bool verbose = false) {
    //     double size_panels = 0;
    //     double size_supp = 0;
    //     for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
    //         size_panels += sdsl::size_in_mega_bytes(phi_vec[i]);
    //         size_panels += sdsl::size_in_mega_bytes(phi_inv_vec[i]);
    //         size_panels += sdsl::size_in_mega_bytes(phi_rank[i]);
    //         size_panels += sdsl::size_in_mega_bytes(phi_select[i]);
    //         size_panels += sdsl::size_in_mega_bytes(phi_inv_rank[i]);
    //         size_panels += sdsl::size_in_mega_bytes(phi_inv_select[i]);
    //         size_supp += sdsl::size_in_mega_bytes(phi_supp[i]);
    //         size_supp += sdsl::size_in_mega_bytes(phi_inv_supp[i]);
    //         size_supp += sdsl::size_in_mega_bytes(phi_supp_l[i]);
    //     }
    //     double size = size_panels + size_supp;
    //     if (verbose) {
    //         std::cout << "phi panels: " << size_panels << " megabytes\n";
    //         std::cout << "phi support: " << size_supp << " megabytes\n";
    //         std::cout << "phi data structure (panels + support): " << size
    //                   << " megabytes\n";
    //     }
    //     return size;
    // }

    /**
     * @brief function to serialize the phi/phi_inv data structure object
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    // size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
    //                  const std::string &name = "") {
    //     sdsl::structure_tree_node *child =
    //             sdsl::structure_tree::add_child(v, name,
    //                                             sdsl::util::class_name(
    //                                                     *this));
    //     size_t written_bytes = 0;
    //     out.write((char *) &this->def, sizeof(this->def));
    //     written_bytes += sizeof(this->def);
    //     out.write((char *) &this->w, sizeof(this->w));
    //     written_bytes += sizeof(this->w);
    //
    //     for (unsigned int i = 0; i < this->phi_vec.size(); i++) {
    //         std::string label = "phi_vec_" + std::to_string(i);
    //         written_bytes += this->phi_vec[i].serialize(out, child, label);
    //     }
    //
    //     for (unsigned int i = 0; i < this->phi_inv_vec.size(); i++) {
    //         std::string label = "phi_inv_vec_" + std::to_string(i);
    //         written_bytes += this->phi_inv_vec[i].serialize(out, child, label);
    //     }
    //
    //     for (unsigned int i = 0; i < this->phi_supp.size(); i++) {
    //         std::string label = "phi_supp_" + std::to_string(i);
    //         written_bytes += this->phi_supp[i].serialize(out, child, label);
    //     }
    //
    //     for (unsigned int i = 0; i < this->phi_inv_supp.size(); i++) {
    //         std::string label = "phi_inv_supp_" + std::to_string(i);
    //         written_bytes += this->phi_inv_supp[i].serialize(out, child,
    //                                                          label);
    //     }
    //
    //     for (unsigned int i = 0; i < this->phi_supp_l.size(); i++) {
    //         std::string label = "phi_supp_l_" + std::to_string(i);
    //         written_bytes += this->phi_supp_l[i].serialize(out, child, label);
    //     }
    //     sdsl::structure_tree::add_size(child, written_bytes);
    //     return written_bytes;
    // }

    /**
     * @brief function to load the phi/phi_inv data structure object
     * @param in std::istream object from which load the phi/phi_inv data
     * structure object
     */
    // void load(std::istream &in) {
    //     in.read((char *) &this->def, sizeof(this->def));
    //     in.read((char *) &this->w, sizeof(this->w));
    //
    //     for (unsigned int i = 0; i < this->def; i++) {
    //         auto s = new sdsl::sd_vector<>();
    //         s->load(in);
    //         this->phi_vec.emplace_back(*s);
    //         delete s;
    //     }
    //     for (unsigned int i = 0; i < this->def; i++) {
    //         auto s = new sdsl::sd_vector<>();
    //         s->load(in);
    //         this->phi_inv_vec.emplace_back(*s);
    //         delete s;
    //     }
    //     for (unsigned int i = 0; i < this->def; i++) {
    //         auto s = new sdsl::int_vector<>();
    //         s->load(in);
    //         this->phi_supp.emplace_back(*s);
    //         delete s;
    //     }
    //     for (unsigned int i = 0; i < this->def; i++) {
    //         auto s = new sdsl::int_vector<>();
    //         s->load(in);
    //         this->phi_inv_supp.emplace_back(*s);
    //         delete s;
    //     }
    //     for (unsigned int i = 0; i < this->def; i++) {
    //         auto s = new sdsl::int_vector<>();
    //         s->load(in);
    //         this->phi_supp_l.emplace_back(*s);
    //         delete s;
    //     }
    //     for (auto &i: this->phi_vec) {
    //         this->phi_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
    //                 &i));
    //         this->phi_select.emplace_back(sdsl::sd_vector<>::select_1_type(
    //                 &i));
    //     }
    //     for (auto &i: this->phi_inv_vec) {
    //         this->phi_inv_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
    //                 &i));
    //         this->phi_inv_select.emplace_back(sdsl::sd_vector<>::select_1_type(
    //                 &i));
    //     }
    // }
};

#endif //PHI_H
