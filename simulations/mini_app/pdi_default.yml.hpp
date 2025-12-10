// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
types:
  tor1_double: &tor1_double
    type: array
    subtype: double
    size: &tor1_size [ '$tor1_extents[0]' ]

  tor2_double: &tor2_double
    type: array
    subtype: double
    size: &tor2_size [ '$tor2_extents[0]' ]

  tor3_double: &tor3_double
    type: array
    subtype: double
    size: &tor3_size [ '$tor3_extents[0]' ]

  vpar_double: &vpar_double
    type: array
    subtype: double
    size: &vpar_size [ '$vpar_extents[0]' ]

  mu_double: &mu_double
    type: array
    subtype: double
    size: &mu_size [ '$mu_extents[0]' ]

metadata:
  gysela_io_filename_size: size_t
  gysela_io_filename: {type: array, subtype: char, size: "$gysela_io_filename_size"}

  #-- Parallel data
  local_fdistribu_starts: { type: array, subtype: size_t, size: 6 }
  local_fdistribu_extents: { type: array, subtype: size_t, size: 6 }

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
  breakpoints_tor3_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_vpar_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_mu_extents: { type: array, subtype: size_t, size: 1 }
  tor1_extents: { type: array, subtype: size_t, size: 1 }
  tor2_extents: { type: array, subtype: size_t, size: 1 }
  tor3_extents: { type: array, subtype: size_t, size: 1 }
  vpar_extents: { type: array, subtype: size_t, size: 1 }
  mu_extents: { type: array, subtype: size_t, size: 1 }
  breakpoints_tor1:
    type: array
    subtype: double
    size: [ '$breakpoints_tor1_extents[0]' ]
  breakpoints_tor2:
    type: array
    subtype: double
    size: [ '$breakpoints_tor2_extents[0]' ]
  breakpoints_tor3:
    type: array
    subtype: double
    size: [ '$breakpoints_tor3_extents[0]' ]
  breakpoints_vpar:
    type: array
    subtype: double
    size: [ '$breakpoints_vpar_extents[0]' ]
  breakpoints_mu:
    type: array
    subtype: double
    size: [ '$breakpoints_mu_extents[0]' ]
  tor1: tor1_double
  tor2: tor2_double
  tor3: tor3_double
  vpar: vpar_double
  mu: mu_double

  # -- Timing statistics
  timing_names_count: size_t
  timing_names_strlen: size_t
  timing_names: { type: array, subtype: char, size: [ '$timing_names_count', '$timing_names_strlen' ] }
  timing_values_size: size_t
  timing_values: { type: array, subtype: double, size: [ '$timing_values_size' ] }

  #-- Fluid moments extents
  density_extents: { type: array, subtype: size_t, size: 4 }
  mean_velocity_extents: { type: array, subtype: size_t, size: 4 }
  temperature_extents: { type: array, subtype: size_t, size: 4 }

data:
  #-- Distribution function
  fdistribu_sptor3Dv2D:
    type: array
    subtype: double
    size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]', '$local_fdistribu_extents[5]' ]

  #-- Fluid moments (4D: species × tor1 × tor2 × tor3)
  density:
    type: array
    subtype: double
    size: [ '$density_extents[0]', '$density_extents[1]', '$density_extents[2]', '$density_extents[3]' ]
  mean_velocity:
    type: array
    subtype: double
    size: [ '$mean_velocity_extents[0]', '$mean_velocity_extents[1]', '$mean_velocity_extents[2]', '$mean_velocity_extents[3]' ]
  temperature:
    type: array
    subtype: double
    size: [ '$temperature_extents[0]', '$temperature_extents[1]', '$temperature_extents[2]', '$temperature_extents[3]' ]

plugins:
  mpi:
  decl_hdf5:
    #-- Read file names and expose them to PDI
    - file: '${gysela_io_filename}'
      on_event: [ReadFileNames]
      read:
        gysela_io_filename_size: ~
        gysela_io_filename: ~
    
    #-- Read species info
    - file: '${gysela_io_filename}'
      on_event: [read_species]
      read:
        species: ~
        masses: ~
        charges: ~

    #-- Read mesh
    - file: '${gysela_io_filename}'
      on_event: [read_tor1_extents]
      read:
        breakpoints_tor1_extents: {size_of: breakpoints_tor1}
        tor1_extents: {size_of: tor1}
    - file: '${gysela_io_filename}'
      on_event: [read_tor1]
      read:
        breakpoints_tor1: ~
        tor1: ~
    - file: '${gysela_io_filename}'
      on_event: [read_tor2_extents]
      read:
        breakpoints_tor2_extents: {size_of: breakpoints_tor2}
        tor2_extents: {size_of: tor2}
    - file: '${gysela_io_filename}'
      on_event: [read_tor2]
      read:
        breakpoints_tor2: ~
        tor2: ~
    - file: '${gysela_io_filename}'
      on_event: [read_tor3_extents]
      read:
        breakpoints_tor3_extents: {size_of: breakpoints_tor3}
        tor3_extents: {size_of: tor3}
    - file: '${gysela_io_filename}'
      on_event: [read_tor3]
      read:
        breakpoints_tor3: ~
        tor3: ~
    - file: '${gysela_io_filename}'
      on_event: [read_vpar_extents]
      read:
        breakpoints_vpar_extents: {size_of: breakpoints_vpar}
        vpar_extents: {size_of: vpar}
    - file: '${gysela_io_filename}'
      on_event: [read_vpar]
      read:
        breakpoints_vpar: ~
        vpar: ~
    - file: '${gysela_io_filename}'
      on_event: [read_mu_extents]
      read:
        breakpoints_mu_extents: {size_of: breakpoints_mu}
        mu_extents: {size_of: mu}
    - file: '${gysela_io_filename}'
      on_event: [read_mu]
      read:
        breakpoints_mu: ~
        mu: ~

    #-- Read distribution function
    - file: '${gysela_io_filename}'
      on_event: [read_fdistribu]
      read:
        fdistribu_sptor3Dv2D:
            dataset_selection:
              size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]', '$local_fdistribu_extents[5]' ]
              start: [ '$local_fdistribu_starts[0]', '$local_fdistribu_starts[1]', '$local_fdistribu_starts[2]', '$local_fdistribu_starts[3]', '$local_fdistribu_starts[4]', '$local_fdistribu_starts[5]' ]

    #-- Write timing (as table: names and values arrays)
    - file: 'cpu_time.h5'
      on_event: [write_timing]
      collision_policy: replace_and_warn
      write:
        timing_names: ~
        timing_values: ~
        
    #-- Write distribution function and coordinates
    - file: 'fdistribu_5D_output.h5'
      on_event: [write_fdistribu]
      collision_policy: replace_and_warn
      communicator: $MPI_COMM_WORLD
      datasets:
        fdistribu_sptor3Dv2D:
          type: array
          subtype: double
          size: [ '$species_extents[0]', '$tor1_extents[0]', '$tor2_extents[0]', '$tor3_extents[0]', '$vpar_extents[0]', '$mu_extents[0]' ]
      write:
        tor1: ~
        tor2: ~
        tor3: ~
        vpar: ~
        mu: ~
        fdistribu_sptor3Dv2D:
            dataset_selection:
              size: [ '$local_fdistribu_extents[0]', '$local_fdistribu_extents[1]', '$local_fdistribu_extents[2]', '$local_fdistribu_extents[3]', '$local_fdistribu_extents[4]', '$local_fdistribu_extents[5]' ]
              start: [ '$local_fdistribu_starts[0]', '$local_fdistribu_starts[1]', '$local_fdistribu_starts[2]', '$local_fdistribu_starts[3]', '$local_fdistribu_starts[4]', '$local_fdistribu_starts[5]' ]

    #-- Write fluid moments
    - file: 'fluid_moments.h5'
      on_event: [write_fluid_moments]
      collision_policy: replace_and_warn
      write:
        density: ~
        mean_velocity: ~
        temperature: ~
trace: ~
)PDI_CFG";
