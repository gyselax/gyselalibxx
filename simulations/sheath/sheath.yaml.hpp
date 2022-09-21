// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PARAMS_CFG(Mesh:
  x_min: 0.0
  x_max: 12.56637061435917
  x_size: 128
  vx_min: -6.0
  vx_max: +6.0
  vx_size: 127

SpeciesInfo:
- charge: -1
  mass: 1.
  density_eq: 1.
  temperature_eq: 1.
  mean_velocity_eq: 0.
  perturb_amplitude: 0.
  perturb_mode: 1

- charge: 1
  mass: 400.
  density_eq: 1.
  temperature_eq: 1.
  mean_velocity_eq: 0.
  perturb_amplitude: 0.
  perturb_mode: 1

SourceInfo:
  px_source: 0.45
  dx_source: 4
  px_sink: 0.45
  dx_sink: 4

KineticSourceInfo:
  source_amplitude: 0.01
  density_amplitude: 1.
  energy_amplitude: 1.
  temperature_source: 1.

Algorithm:
  deltat: 0.125
  nbiter: 360

Output:
  time_diag: 0.25
)PARAMS_CFG";
