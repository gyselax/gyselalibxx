// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nvpar_spline_cells : int
  Nmu_spline_cells : int
  iter : int
  nbstep_diag: int
  iter_start : int
  iter_saved : int
  time_saved : double
  grid_vpar_extents: { type: array, subtype: int64, size: 1 }
  grid_vpar:
    type: array
    subtype: double
    size: [ '$grid_vpar_extents[0]' ]
  grid_mu_extents: { type: array, subtype: int64, size: 1 }
  grid_mu:
    type: array
    subtype: double
    size: [ '$grid_mu_extents[0]' ]
  fdistribu_charges_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_charges:
    type: array
    subtype: double
    size: [ '$fdistribu_charges_extents[0]' ]
  fdistribu_masses_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_masses:
    type: array
    subtype: double
    size: [ '$fdistribu_masses_extents[0]' ]
  fdistribu_eq_extents : { type: array, subtype: int64, size: 3 }
  fdistribu_eq:
    type: array
    subtype: double
    size: [ '$fdistribu_eq_extents[0]', '$fdistribu_eq_extents[1]', '$fdistribu_eq_extents[2]' ]

data:
  fdistribu_extents: { type: array, subtype: int64, size: 3 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]', '$fdistribu_extents[2]' ]

plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iter:
        - set:
          - iter_saved: '${iter_start} + ${iter}/${nbstep_diag}'
    on_finalise:
      - release: [iter_saved]

  decl_hdf5:
    - file: 'coll_${iter_saved:05}.h5'
      on_event: [write_fdistribu]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [Nvpar_spline_cells, Nmu_spline_cells, grid_vpar, grid_mu, nbstep_diag, fdistribu_charges, fdistribu_masses, time_saved, fdistribu]
  #trace: ~
)PDI_CFG";
