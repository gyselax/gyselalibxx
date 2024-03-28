// SPDX-License-Identifier: MIT

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  read_restart_filename_size: size_t
  read_restart_filename: {type: array, subtype: char, size: "$read_restart_filename_size"}
  write_restart_filename_size: size_t
  write_restart_filename: {type: array, subtype: char, size: "$write_restart_filename_size"}
  iter_start : int
  iter_saved : int
  time_saved : double
  grid_tor1_extents: { type: array, subtype: size_t, size: 1 }
  grid_tor2_extents: { type: array, subtype: size_t, size: 1 }
  grid_tor3_extents: { type: array, subtype: size_t, size: 1 }
  grid_vpar_extents: { type: array, subtype: size_t, size: 1 }
  grid_mu_extents: { type: array, subtype: size_t, size: 1 }
  species_extents: { type: array, subtype: size_t, size: 1 }
  masses_extents: { type: array, subtype: size_t, size: 1 }
  charges_extents: { type: array, subtype: size_t, size: 1 }
data:
  grid_tor1:
    type: array
    subtype: double
    size: [ '$grid_tor1_extents[0]' ]
  grid_tor2:
    type: array
    subtype: double
    size: [ '$grid_tor2_extents[0]' ]
  grid_tor3:
    type: array
    subtype: double
    size: [ '$grid_tor3_extents[0]' ]
  grid_vpar:
    type: array
    subtype: double
    size: [ '$grid_vpar_extents[0]' ]
  grid_mu:
    type: array
    subtype: double
    size: [ '$grid_mu_extents[0]' ]
  species:
    type: array
    subtype: int
    size: [ '$species_extents[0]' ]
  masses:
    type: array
    subtype: double
    size: [ '$masses_extents[0]' ]
  charges:
    type: array
    subtype: double
    size: [ '$charges_extents[0]' ]
  densityTorCS:
    type: array
    subtype: double
    size: [ '$species_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]
  temperatureTorCS:
    type: array
    subtype: double
    size: [ '$species_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]
  UparTorCS:
    type: array
    subtype: double
    size: [ '$species_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]
  fdistribu:
    type: array
    subtype: double
    size: [ '$species_extents[0]', '$grid_tor3_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]', '$grid_vpar_extents[0]', '$grid_mu_extents[0]' ]
 
plugins:
  decl_hdf5:
    - file: '${read_restart_filename}'
      on_event: [read_grid_extents]
      read:
        grid_tor1_extents: {size_of: grid_tor1}
        grid_tor2_extents: {size_of: grid_tor2}
        grid_tor3_extents: {size_of: grid_tor3}
        grid_vpar_extents: {size_of: grid_vpar}
        grid_mu_extents: {size_of: grid_mu}
        species_extents: {size_of: species}
        masses_extents: {size_of: masses}
        charges_extents: {size_of: charges}
    - file: '${read_restart_filename}'
      on_event: [read_grid]
      read:
        grid_tor1: ~
        grid_tor2: ~
        grid_tor3: ~
        grid_vpar: ~
        grid_mu: ~
        species: ~
        masses: ~
        charges: ~
    - file: '${read_restart_filename}'
      on_event: [read_profiles]
      read:
        - densityTorCS
        - temperatureTorCS
        - UparTorCS
    - file: '${read_restart_filename}'
      on_event: [read_fdistribu]
      read:
        - time_saved
        - fdistribu
    - file: '${write_restart_filename}'
      on_event: [write_restart]
      collision_policy: skip_and_warn
      write:
       - iter_saved
       - time_saved
       - grid_tor1
       - grid_tor2
       - grid_tor3
       - grid_vpar
       - grid_mu
       - species
       - masses
       - charges
       - densityTorCS
       - temperatureTorCS
       - UparTorCS
       - fdistribu
trace: ~
)PDI_CFG";
