// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  r_ncells: 128
  r_min: 0.0
  r_max: 1.0 
  p_ncells: 256
  p_min: 0.0
  p_max: 6.283185307179586

Time:
  delta_t: 0.1
  final_T: 10.0
  
Perturbation:
  eps: 0.0001
  x_star_1: 0.08
  y_star_1: -0.14
  x_star_2: -0.08
  y_star_2: 0.14
  sigma: 0.08
  
Output:
  time_step_diag: 5
)PDI_CFG";
