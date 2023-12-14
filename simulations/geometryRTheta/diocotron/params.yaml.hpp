// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(Mesh:
  r_size: 128
  p_size: 256
  r_min: 0.0
  r_minus: 0.45
  r_plus: 0.50
  r_max: 1.0 

Time:
  delta_t: 0.1
  final_T: 100.0
  
Perturbation:
  charge_Q: 0
  l_mode: 9 
  eps: 0.0001
  
Output:
  time_step_diag: 10
)PDI_CFG";
