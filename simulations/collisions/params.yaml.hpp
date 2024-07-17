// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(SplineMesh:
  vpar_min: -6.0
  vpar_max: +6.0
  vpar_ncells: 127
  mu_min: 0.0
  mu_max: 12.0
  mu_ncells: 63

SpeciesInfo:
- charge: 1
  mass: 1.
  density_eq: 1.
  temperature_eq: 1.
  mean_velocity_eq: 0.
- charge: -1
  mass: 5.44e-4
  density_eq: 1.
  temperature_eq: 1.1
  mean_velocity_eq: 0.

Algorithm:
  deltat: 0.01
  nbiter: 10

Output:
  time_diag: 0.1
)PDI_CFG";
