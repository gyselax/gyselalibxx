// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  r_size: int
  p_size: int
  r_min: double
  r_max: double 

  iter_saved: int
  time_step_diag: int

  iter: int
  time: double
  delta_t: double
  final_T: double

  eps: double
  sigma: double
  x_star_1: double 
  y_star_1: double
  x_star_2: double
  y_star_2: double


  r_coords_extents: {type: array,  subtype: int64, size: 1 }
  r_coords:
    type: array
    subtype: double
    size: [  '$r_coords_extents[0]' ]

  p_coords_extents: {type: array,  subtype: int64, size: 1 }
  p_coords:
    type: array
    subtype: double
    size: [  '$p_coords_extents[0]' ]




data:
  x_coords_extents: {type: array,  subtype: int64, size: 2 }
  x_coords:
    type: array
    subtype: double
    size: [  '$x_coords_extents[0]', '$x_coords_extents[1]' ]

  y_coords_extents: {type: array,  subtype: int64, size: 2 }
  y_coords:
    type: array
    subtype: double
    size: [  '$y_coords_extents[0]', '$y_coords_extents[1]' ]



  jacobian_extents: {type: array,  subtype: int64, size: 2 }
  jacobian: 
    type: array
    subtype: double
    size: [  '$jacobian_extents[0]', '$jacobian_extents[1]' ]


  density_eq_extents: {type: array, subtype: int64, size: 2 }
  density_eq: 
    type: array
    subtype: double
    size: [ '$density_eq_extents[0]', '$density_eq_extents[1]' ]

  electrical_potential_eq_extents: {type: array, subtype: int64, size: 2 }
  electrical_potential_eq: 
    type: array
    subtype: double
    size: [ '$electrical_potential_eq_extents[0]', '$electrical_potential_eq_extents[1]' ]

  density_extents: {type: array, subtype: int64, size: 2 }
  density: 
    type: array
    subtype: double
    size: [ '$density_extents[0]', '$density_extents[1]' ]


  electrical_potential_extents: {type: array, subtype: int64, size: 2 }
  electrical_potential: 
    type: array
    subtype: double
    size: [ '$electrical_potential_extents[0]', '$electrical_potential_extents[1]' ]



plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iteration:
        - set:
          - iter_saved: '${iter} % ${time_step_diag}' 
    on_finalize:
      - release: [iter_saved]

  decl_hdf5:
    - file: 'output/VOICEXX_initstate.h5'
      on_event: [initialization]
      collision_policy: replace_and_warn
      write: [r_size, p_size, r_coords, p_coords, x_coords, y_coords, jacobian, delta_t, final_T, time_step_diag, density_eq, electrical_potential_eq]

    - file: 'output/VOICEXX_${iter:05}.h5'
      on_event: [iteration, last_iteration]
      when: '${iter} % ${time_step_diag} = 0'
      collision_policy: replace_and_warn
      
      write: [time, density, electrical_potential]
  #trace: ~
)PDI_CFG";
