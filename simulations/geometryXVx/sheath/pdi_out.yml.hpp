// SPDX-License-Identifier: MIT

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx : int
  Nvx : int
  iter : int
  iter_start : int
  time_saved : double
  nbstep_diag: int
  iter_saved : int
  Lx : double
  MeshX_extents: { type: array, subtype: int64, size: 1 }
  MeshX:
    type: array
    subtype: double
    size: [ '$MeshX_extents[0]' ]
  MeshVx_extents: { type: array, subtype: int64, size: 1 }
  MeshVx:
    type: array
    subtype: double
    size: [ '$MeshVx_extents[0]' ]
  Nkinspecies: int
  fdistribu_charges_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_charges:
    type: array
    subtype: int
    size: [ '$fdistribu_charges_extents[0]' ]
  fdistribu_masses_extents : { type: array, subtype: int64, size: 1 }
  fdistribu_masses:
    type: array
    subtype: double
    size: [ '$fdistribu_masses_extents[0]' ]
  fdistribu_eq_extents : { type: array, subtype: int64, size: 2 }
  fdistribu_eq:
    type: array
    subtype: double
    size: [ '$fdistribu_eq_extents[0]', '$fdistribu_eq_extents[1]' ]
  collintra_nustar0 : double
  collinter_nustar0 : double

  krook_sink_adaptive_extent : double
  krook_sink_adaptive_stiffness : double
  krook_sink_adaptive_amplitude : double
  krook_sink_adaptive_density : double
  krook_sink_adaptive_temperature : double
  krook_sink_adaptive_mask_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_adaptive_mask:
    type: array
    subtype: double
    size: [ '$krook_sink_adaptive_mask_extents[0]' ]
  krook_sink_adaptive_ftarget_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_adaptive_ftarget:
    type: array
    subtype: double
    size: [ '$krook_sink_adaptive_ftarget_extents[0]' ]
  
  krook_sink_constant_extent : double
  krook_sink_constant_stiffness : double
  krook_sink_constant_amplitude : double
  krook_sink_constant_density : double
  krook_sink_constant_temperature : double
  krook_sink_constant_mask_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_constant_mask:
    type: array
    subtype: double
    size: [ '$krook_sink_constant_mask_extents[0]' ]
  krook_sink_constant_ftarget_extents: { type: array, subtype: int64, size: 1 }
  krook_sink_constant_ftarget:
    type: array
    subtype: double
    size: [ '$krook_sink_constant_ftarget_extents[0]' ]

  kinetic_source_extent : double
  kinetic_source_stiffness : double
  kinetic_source_amplitude : double
  kinetic_source_density : double
  kinetic_source_energy : double
  kinetic_source_temperature : double
  kinetic_source_velocity_shape_extents: { type: array, subtype: int64, size: 1 }
  kinetic_source_velocity_shape:
    type: array
    subtype: double
    size: [ '$kinetic_source_velocity_shape_extents[0]' ]
  kinetic_source_spatial_extent_extents: { type: array, subtype: int64, size: 1 }
  kinetic_source_spatial_extent:
    type: array
    subtype: double
    size: [ '$kinetic_source_spatial_extent_extents[0]' ]

data:
  fdistribu_extents: { type: array, subtype: int64, size: 3 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]', '$fdistribu_extents[2]' ]
  electrostatic_potential_extents: { type: array, subtype: int64, size: 1 }
  electrostatic_potential:
    type: array
    subtype: double
    size: [ '$electrostatic_potential_extents[0]' ]

plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iter:
        - set:
          - iter_saved: '${iter_start} + ${iter}/${nbstep_diag}'
    on_finalize:
      - release: [iter_saved]
  decl_hdf5:
    - file: 'VOICEXX_initstate.h5'
      on_event: [initial_state]
      collision_policy: replace_and_warn
      write:
        - Nx
        - Nvx
        - MeshX
        - MeshVx
        - Lx
        - nbstep_diag
        - collintra_nustar0
        - collinter_nustar0
        - Nkinspecies
        - fdistribu_charges
        - fdistribu_masses
        - fdistribu_eq

        - krook_sink_adaptive_extent
        - krook_sink_adaptive_stiffness
        - krook_sink_adaptive_amplitude
        - krook_sink_adaptive_density
        - krook_sink_adaptive_temperature
        - krook_sink_adaptive_mask
        - krook_sink_adaptive_ftarget

        - krook_sink_constant_extent
        - krook_sink_constant_stiffness
        - krook_sink_constant_amplitude
        - krook_sink_constant_density
        - krook_sink_constant_temperature
        - krook_sink_constant_mask
        - krook_sink_constant_ftarget

        - kinetic_source_extent
        - kinetic_source_stiffness
        - kinetic_source_amplitude
        - kinetic_source_density
        - kinetic_source_energy
        - kinetic_source_temperature
        - kinetic_source_velocity_shape
        - kinetic_source_spatial_extent
    - file: 'VOICEXX_${iter_saved:05}.h5'
      on_event: [iteration, last_iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [time_saved, fdistribu, electrostatic_potential]
    - file: 'VOICEXX_${iter_start:05}.h5'
      on_event: restart
      read: [time_saved, fdistribu]
  #trace: ~
)PDI_CFG";
