// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  read_restart_filename_size: size_t
  read_restart_filename: {type: array, subtype: char, size: "$read_restart_filename_size"}
  write_restart_filename_size: size_t
  write_restart_filename: {type: array, subtype: char, size: "$write_restart_filename_size"}
  iter_start : int
  iter_saved : int
  time_saved : double

  #-- Parallel data
  local_fdistribu_starts: { type: array, subtype: size_t, size: 5 }
  local_fdistribu_extents: { type: array, subtype: size_t, size: 5 }

  #-- Species info
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
  species_extents: { type: array, subtype: size_t, size: 1 }
  masses_extents: { type: array, subtype: size_t, size: 1 }
  charges_extents: { type: array, subtype: size_t, size: 1 }

  #-- Mesh: breakpoints and grid
  breakpoints_tor1_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_tor2_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_vpar_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_mu_extents: { type: array, subtype: size_t, size: 1 }
  grid_tor1_extents: { type: array, subtype: size_t, size: 1 }
  grid_tor2_extents: { type: array, subtype: size_t, size: 1 }
  grid_vpar_extents: { type: array, subtype: size_t, size: 1 }
  grid_mu_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_tor1:
    type: array
    subtype: double
    size: [ '$breakpoints_tor1_extents[0]' ]
  breakpoints_tor2:
    type: array
    subtype: double
    size: [ '$breakpoints_tor2_extents[0]' ]
  breakpoints_vpar:
    type: array
    subtype: double
    size: [ '$breakpoints_vpar_extents[0]' ]
  breakpoints_mu:
    type: array
    subtype: double
    size: [ '$breakpoints_mu_extents[0]' ]
  grid_tor1:
    type: array
    subtype: double
    size: [ '$grid_tor1_extents[0]' ]
  grid_tor2:
    type: array
    subtype: double
    size: [ '$grid_tor2_extents[0]' ]
  grid_vpar:
    type: array
    subtype: double
    size: [ '$grid_vpar_extents[0]' ]
  grid_mu:
    type: array
    subtype: double
    size: [ '$grid_mu_extents[0]' ]

  #-- Magnetic config: R_matrix and Z_matrix
  R_matrix:
    type: array
    subtype: double
    size: [ '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]
  Z_matrix:
    type: array
    subtype: double
    size: [ '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]
  
  #-- Magnetic config: normB_matrix
  normB_matrix:
    type: array
    subtype: double
    size: [ '$grid_tor2_extents[0]', '$grid_tor1_extents[0]' ]

data:
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
    size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]' ]
 
plugins:
  mpi:
  decl_hdf5:
    #-- Read species info
    - file: '${read_restart_filename}'
      on_event: [read_species_extents]
      read:
        species_extents: {size_of: species}
        masses_extents: {size_of: masses}
        charges_extents: {size_of: charges}
    - file: '${read_restart_filename}'
      on_event: [read_species]
      read:
        species: ~
        masses: ~
        charges: ~
    
    #-- Read mesh
    - file: '${read_restart_filename}'
      on_event: [read_tor1_extents]
      read:
        breakpoints_tor1_extents: {size_of: breakpoints_tor1}
        grid_tor1_extents: {size_of: grid_tor1}
    - file: '${read_restart_filename}'
      on_event: [read_tor1]
      read:
        breakpoints_tor1: ~
        grid_tor1: ~
    - file: '${read_restart_filename}'
      on_event: [read_tor2_extents]
      read:
        breakpoints_tor2_extents: {size_of: breakpoints_tor2}
        grid_tor2_extents: {size_of: grid_tor2}
    - file: '${read_restart_filename}'
      on_event: [read_tor2]
      read:
        breakpoints_tor2: ~
        grid_tor2: ~
    - file: '${read_restart_filename}'
      on_event: [read_vpar_extents]
      read:
        breakpoints_vpar_extents: {size_of: breakpoints_vpar}
        grid_vpar_extents: {size_of: grid_vpar}
    - file: '${read_restart_filename}'
      on_event: [read_vpar]
      read:
        breakpoints_vpar: ~
        grid_vpar: ~
    - file: '${read_restart_filename}'
      on_event: [read_mu_extents]
      read:
        breakpoints_mu_extents: {size_of: breakpoints_mu}
        grid_mu_extents: {size_of: grid_mu}
    - file: '${read_restart_filename}'
      on_event: [read_mu]
      read:
        breakpoints_mu: ~
        grid_mu: ~

    #-- Read magnetic config
    - file: '${read_restart_filename}'
      on_event: [read_magnetic_config]
      read:
        R_matrix: ~
        Z_matrix: ~
        normB_matrix: ~
    
    #-- Read profiles
    - file: '${read_restart_filename}'
      on_event: [read_profiles]
      read:
        - densityTorCS
        - temperatureTorCS
        - UparTorCS

    #-- Read full distribution function
    - file: '${read_restart_filename}'
      on_event: [read_fdistribu]
      communicator: $MPI_COMM_WORLD
      datasets:
        fdistribu:
          type: array
          subtype: double
          size: [ '$species_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]', '$grid_vpar_extents[0]', '$grid_mu_extents[0]' ]
      read:
        time_saved: ~
        fdistribu:
          dataset_selection:
            size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]' ]
            start: [ '$local_fdistribu_starts[0]', '$local_fdistribu_starts[1]', '$local_fdistribu_starts[2]', '$local_fdistribu_starts[3]', '$local_fdistribu_starts[4]' ]

    - file: '${write_restart_filename}'
      on_event: [write_restart]
      collision_policy: skip_and_warn
      communicator: $MPI_COMM_WORLD
      datasets:
        fdistribu:
          type: array
          subtype: double
          size: [ '$species_extents[0]', '$grid_tor2_extents[0]', '$grid_tor1_extents[0]', '$grid_vpar_extents[0]', '$grid_mu_extents[0]' ]
      write:
       iter_saved: ~
       time_saved: ~
       species: ~
       masses: ~
       charges: ~
       breakpoints_tor1: ~
       breakpoints_tor2: ~
       breakpoints_vpar: ~
       breakpoints_mu: ~
       grid_tor1: ~
       grid_tor2: ~
       grid_vpar: ~
       grid_mu: ~
       R_matrix: ~
       Z_matrix: ~
       normB_matrix: ~
       densityTorCS: ~
       temperatureTorCS: ~
       UparTorCS: ~
       fdistribu:
         dataset_selection:
           size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]' ]
           start: [ '$local_fdistribu_starts[0]', '$local_fdistribu_starts[1]', '$local_fdistribu_starts[2]', '$local_fdistribu_starts[3]', '$local_fdistribu_starts[4]' ]
trace: ~
)PDI_CFG";
