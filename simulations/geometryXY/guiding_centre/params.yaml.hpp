// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  x_min: 0.0
  x_max: 12.56637061435917
  x_ncells: 64
  y_min: 0.0
  y_max: 12.56637061435917
  y_ncells : 64

PerturbationInfo:
  perturb_amplitude: 0.015
  perturb_mode: 0.5

Algorithm:
  delta_t: 0.05
  final_time: 30

Output:
  nbstep_diag: 4
)PDI_CFG";
