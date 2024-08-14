// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  r_min: 0.
  r_max: 1.
  r_ncells: 128
  p_ncells: 256
Time:
  time_step: 1e-2
  final_time: 1
Output:
  save_curves: false
  save_feet: false
)PDI_CFG";
