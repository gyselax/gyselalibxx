constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx : int
  Nvx : int
  iter : int
  iter_time : double
  nbstep_diag: int
  iter_save : int
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
  
data:
  fdistribu_extents: { type: array, subtype: int64, size: 3 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]', '$fdistribu_extents[2]' ]
  efield_extents: { type: array, subtype: int64, size: 1 }
  efield:
    type: array
    subtype: double
    size: [ '$efield_extents[0]' ]
    
plugins:
  set_value:
    on_init:
      - share:
        - iter_save: 0
    on_data:
      iter:
        - set:
          - iter_save: '${iter}/${nbstep_diag}'
    on_finalize:
      - release: [iter_save]
  decl_hdf5:    
    - file: 'VOICEXX_initstate.h5'
      on_event: [initial_state]
      collision_policy: replace_and_warn
      write: [Nx, Nvx, MeshX, MeshVx, nbstep_diag ]
    - file: 'VOICEXX_${iter_save:05}.h5'
      on_event: [iteration, last_iteration]
      when: '${iter} % ${nbstep_diag} = 0'
      collision_policy: replace_and_warn
      write: [iter, iter_time, Nx, Nvx, MeshX, MeshVx, fdistribu, efield ]
  #trace: ~
)PDI_CFG";
