// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx_spline_cells : int
  Ny_spline_cells : int
  Nvx_spline_cells : int
  Nvy_spline_cells : int
  iter : int
  time_saved : double
  nbstep_diag: int
  iter_saved : int
  MeshX_extents: { type: array, subtype: int64, size: 1 }
  MeshX:
    type: array
    subtype: double
    size: [ '$MeshX_extents[0]' ]
  MeshY_extents: { type: array, subtype: int64, size: 1 }
  MeshY:
    type: array
    subtype: double
    size: [ '$MeshY_extents[0]' ]
  MeshVx_extents: { type: array, subtype: int64, size: 1 }
  MeshVx:
    type: array
    subtype: double
    size: [ '$MeshVx_extents[0]' ]
  MeshVy_extents: { type: array, subtype: int64, size: 1 }
  MeshVy:
    type: array
    subtype: double
    size: [ '$MeshVy_extents[0]' ]
  Nkinspecies: int
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

  #-- Parallel data
  local_fdistribu_starts: { type: array, subtype: size_t, size: 5 }
  local_fdistribu_extents: { type: array, subtype: size_t, size: 5 }


data:
  fdistribu:
    type: array
    subtype: double
    size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]' ]
  electrostatic_potential_extents: { type: array, subtype: int64, size: 2 }
  electrostatic_potential:
    type: array
    subtype: double
    size: [ '$electrostatic_potential_extents[0]', '$electrostatic_potential_extents[1]' ]

plugins:
  mpi:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iter:
        - set:
          - iter_saved: '${iter}/${nbstep_diag}'
    on_finalize:
      - release: [iter_saved]
  decl_hdf5:
    - file: 'GYSELALIBXX_initstate.h5'
      on_event: [initial_state]
      collision_policy: replace_and_warn
      write: [Nx_spline_cells, Nvx_spline_cells, MeshX, MeshY, MeshVx, MeshVy, nbstep_diag, Nkinspecies, fdistribu_charges, fdistribu_masses, fdistribu_eq]
    - file: 'GYSELALIBXX_${iter_saved:05}.h5'
      communicator: $MPI_COMM_WORLD
      on_event: [iteration, last_iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      datasets:
        fdistribu:
          type: array
          subtype: double
          size: [ '$Nkinspecies', '$MeshX_extents[0]', '$MeshY_extents[0]', '$MeshVx_extents[0]', '$MeshVy_extents[0]' ]
      write:
        time_saved: ~
        fdistribu:
          dataset_selection:
            size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]' ]
            start: [ '$local_fdistribu_starts[0]', '$local_fdistribu_starts[1]', '$local_fdistribu_starts[2]', '$local_fdistribu_starts[3]', '$local_fdistribu_starts[4]' ]
        electrostatic_potential: ~
  #trace: ~
)PDI_CFG";
