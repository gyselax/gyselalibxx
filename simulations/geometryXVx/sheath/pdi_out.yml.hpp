constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx : int
  Nvx : int
  iter : int
  time_saved : double
  nbstep_diag: int
  iter_saved : int
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

  collintra_nustar0 : double

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
          - iter_saved: '${iter}/${nbstep_diag}'
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
        - nbstep_diag
        - collintra_nustar0
    - file: 'VOICEXX_${iter_saved:05}.h5'
      on_event: [iteration, last_iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [time_saved, fdistribu, electrostatic_potential ]
  #trace: ~
)PDI_CFG";
