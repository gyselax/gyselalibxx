// SPDX-License-Identifier: MIT

#pragma once

constexpr char const* const params_yaml = R"PDI_CFG(
InputFileNames:
  read_restart : GyselaX_restart5D_00000.h5
  write_restart : GyselaX_restart5D_00001.h5
 
Algorithm:
  deltat: 0.125
  nbiter: 360

Output:
  time_diag: 0.25
)PDI_CFG";
