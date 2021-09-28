constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  Nx : int
  Nvx : int
  iter : int
  iter_time : double
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
  fdistribu_extents: { type: array, subtype: int64, size: 2 }
  fdistribu:
    type: array
    subtype: double
    size: [ '$fdistribu_extents[0]', '$fdistribu_extents[1]' ]
  ex_extents: { type: array, subtype: int64, size: 1 }
  ex:
    type: array
    subtype: double
    size: [ '$ex_extents[0]' ]
    
plugins:
  decl_hdf5:
    - file: 'VOICEXX_${iter:05}.h5'
      on_event: [iteration, last_iteration]
      collision_policy: replace_and_warn
      write: [iter, iter_time, Nx, Nvx, MeshX, MeshVx, fdistribu, ex ]
  trace: ~
)PDI_CFG";
