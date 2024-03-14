// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  x_min: 0.0
  x_max: 50.
  x_ncells: 256
  vx_min: -8.0
  vx_max: +8.0
  vx_ncells: 127

SpeciesInfo:
- charge: -1
  mass: 0.0005
  epsilon_bot: 0.1
  temperature_bot: 0.2 # [ .1 ; .5 ]
  mean_velocity_bot: 3.8 # [ 2 ; 6 ]
  perturb_amplitude: 0.00001 # [ 0 ; 0.001 ]
  perturb_mode: 3

Algorithm:
  deltat: 0.1
  nbiter: 600

Output:
  time_diag: 0.4
)PDI_CFG";

// Question: BOT with several species ?
