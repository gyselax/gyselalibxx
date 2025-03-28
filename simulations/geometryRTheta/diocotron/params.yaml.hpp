// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  r_ncells: 128
  r_min: 0.0
  r_max: 1.0 
  theta_ncells: 256
  theta_min: 0.0
  theta_max: 6.283185307179586

Time:
  delta_t: 0.1
  final_T: 100.
  
Perturbation:
  charge_Q: 0.
  l_mode: 9 
  eps: 0.0001
  r_min: 0.45
  r_max: 0.50
  
Output:
  time_step_diag: 10
)PDI_CFG";
