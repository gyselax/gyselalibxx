// SPDX-License-Identifier: MIT

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx_spline_cells : int
  Ny_spline_cells : int

  iter : int
  nbstep_diag: int
  time_step : double
  time_saved : double
  final_time : double

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

  fdistribu_equilibrium_extents : { type: array, subtype: int64, size: 2 }
  fdistribu_equilibrium:
    type: array
    subtype: double
    size: [ '$fdistribu_equilibrium_extents[0]', '$fdistribu_equilibrium_extents[1]' ]

data:
  fdistribu_extents : { type: array, subtype: int64, size: 2 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]' ]


  electrostatic_potential_extents: { type: array, subtype: int64, size: 2 }
  electrostatic_potential:
    type: array
    subtype: double
    size: [ '$electrostatic_potential_extents[0]', '$electrostatic_potential_extents[1]' ]


  electric_field_x_extents: { type: array, subtype: int64, size: 2 }
  electric_field_x:
    type: array
    subtype: double
    size: [ '$electric_field_x_extents[0]', '$electric_field_x_extents[1]' ]

  electric_field_y_extents: { type: array, subtype: int64, size: 2 }
  electric_field_y:
    type: array
    subtype: double
    size: [ '$electric_field_y_extents[0]', '$electric_field_y_extents[1]' ]


plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iteration:
        - set:
          - iter_saved:  '${iter} % ${nbstep_diag}' 
    on_finalize:
      - release: [iter_saved]
      
  decl_hdf5:
    - file: 'output/VOICEXX_initstate.h5'
      on_event: [initialization]
      collision_policy: replace_and_warn
      write: [Nx_spline_cells, Ny_spline_cells, MeshX, MeshY, nbstep_diag, final_time, time_step, fdistribu_equilibrium]

    - file: 'output/VOICEXX_${iter:05}.h5'
      on_event: [iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [iter, time_saved, fdistribu, electrostatic_potential, electric_field_x, electric_field_y]
  #trace: ~
)PDI_CFG";
