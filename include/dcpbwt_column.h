//
// Created by shakyap on 8/12/24.
//

#ifndef DCPBWT_COLUMN_H
#define DCPBWT_COLUMN_H

#include <iostream>
#include "dynamic/dynamic.hpp"

using namespace dyn;
using namespace std;

class dcpbwt_column {
 public:
  packed_spsi zeros;
  packed_spsi ones;
  packed_spsi combined;
//  packed_spsi pref_samples_beg;
//  packed_spsi pref_samples_end;
  packed_spsi pref_samples_beg;
  packed_spsi pref_samples_end;
  packed_spsi div_samples_beg;
  bool start_with_zero;
  // TODO: don't need this, as zeros can get this number
  // TODO: For now keeping it for the sake of ease
  unsigned int num_zeros;

  dcpbwt_column(packed_spsi zeros,
                packed_spsi ones,
                packed_spsi combined,
                packed_spsi pref_samples_beg,
                packed_spsi pref_samples_end,
                packed_spsi div_samples_beg,
                bool start_with_zero,
                unsigned int num_zeros) : zeros(std::move(zeros)),
                                          ones(std::move(ones)),
                                          combined(std::move(combined)),
                                          pref_samples_beg(std::move(pref_samples_beg)),
                                          pref_samples_end(std::move(pref_samples_end)),
                                          div_samples_beg(std::move(div_samples_beg)),
                                          start_with_zero(start_with_zero),
                                          num_zeros(num_zeros) {}

};

#endif //DCPBWT_COLUMN_H
