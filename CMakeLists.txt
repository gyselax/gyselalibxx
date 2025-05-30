# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.20)

# Read version from external file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/src/version.txt" GYSELALIBXX_VERSION_STRING)
string(STRIP "${GYSELALIBXX_VERSION_STRING}" GYSELALIBXX_VERSION_STRING)

# Split into components
string(REPLACE "." ";" VERSION_LIST ${GYSELALIBXX_VERSION_STRING})
list(GET VERSION_LIST 0 GYSELALIBXX_VERSION_MAJOR)
list(GET VERSION_LIST 1 GYSELALIBXX_VERSION_MINOR)
list(GET VERSION_LIST 2 GYSELALIBXX_VERSION_PATCH)

# Set project language and version
project(gyselalibxx
    LANGUAGES C CXX
    VERSION ${GYSELALIBXX_VERSION_MAJOR}.${GYSELALIBXX_VERSION_MINOR}.${GYSELALIBXX_VERSION_PATCH})

set_property(GLOBAL PROPERTY GLOBAL_DEPENDS_NO_CYCLES ON)

# Configuration

## We use C++17 all around!
set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "The C++ standard whose features are requested to build this project.")

## List of options
option(GYSELALIBXX_BUILD_SIMULATIONS "Build the simulations." ON)
option(GYSELALIBXX_BUILD_TESTING "Build the tests." ON)
option(GYSELALIBXX_ENABLE_DEPRECATED "Enable deprecated code" OFF)
option(GYSELALIBXX_COMPILE_SOURCE "Enable compilation of the source code (this should be set to off to build documentation without the C++ dependencies)" ON)
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "-O1" CACHE STRING "Default flags for C++ specific to Gyselalib++")
option(ACTIVATE_RESTART_TESTS "Activate tests which check that a simulation gives the same results after restart." ON)

set(GYSELALIBXX_DEPENDENCY_POLICIES "AUTO" "EMBEDDED" "INSTALLED")

set(GYSELALIBXX_INCLUDE_TESTING_DEPENDENCIES ${GYSELALIBXX_BUILD_TESTING} CACHE BOOL "Indicates if the testing dependencies should be included. By default they are only included if the tests are built.")

# Set default DDC options when included in GYSELA
option(DDC_BUILD_BENCHMARKS       "Build DDC benchmarks." OFF)
option(DDC_BUILD_EXAMPLES         "Build DDC examples" OFF)
option(DDC_BUILD_TESTS            "Build DDC tests if BUILD_TESTING is enabled" OFF)

###############################################################################################
#                                     Get dependencies
###############################################################################################

if ("${GYSELALIBXX_COMPILE_SOURCE}")

# Our dependencies

## Use CTest for running tests
include(CTest)

## Look for a pre-installed LAPACK
find_package(LAPACK REQUIRED)

## Look for a pre-installed MPI
find_package(MPI REQUIRED)

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)

## Look for a pre-installed PDI
find_package(PDI REQUIRED COMPONENTS C)

## if tests are enabled, use googletest from `vendor/`
### we use it to write unit tests
if("${GYSELALIBXX_INCLUDE_TESTING_DEPENDENCIES}")
  set(GYSELALIBXX_GTest_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `GTest` package. Options: ${GYSELALIBXX_DEPENDENCY_POLICIES}")
  set_property(CACHE GYSELALIBXX_GTest_DEPENDENCY_POLICY PROPERTY STRINGS ${GYSELALIBXX_DEPENDENCY_POLICIES})

  if("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "AUTO")
    find_package(GTest "1.12" QUIET)
    if(NOT "${GTest_FOUND}")
      add_subdirectory("vendor/googletest/" "googletest") # SYSTEM)
    endif()
  elseif("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
    add_subdirectory("vendor/googletest/" "googletest") # SYSTEM)
  elseif("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
    find_package(GTest "1.12" REQUIRED)
  endif()
endif()

find_package(Kokkos 4.5.1 QUIET)
if(NOT "${Kokkos_FOUND}")
  ## Use Kokkos from `vendor/`
  if("${Kokkos_ENABLE_CUDA}")
    option(Kokkos_ENABLE_CUDA_CONSTEXPR "Whether to activate experimental relaxed constexpr functions" ON)
    option(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for CUDA" ON)
  endif()
  add_subdirectory("vendor/kokkos/" "kokkos") # SYSTEM)
endif()

find_package(KokkosKernels 4.5.1 QUIET)
if(NOT "${KokkosKernels_FOUND}")
  # NOTE: Do not forget to keep Kokkos and Kokkos-kernels in sync !
  set(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "")
  set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "")
  set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "")
  add_subdirectory("vendor/kokkos-kernels/" "kokkos-kernels") # SYSTEM)
endif()


find_package(DDC 0.5 QUIET COMPONENTS fft pdi splines)
if(NOT "${DDC_FOUND}")
  ## Use the discrete domain computation library (ddc) from `vendor/`
  add_subdirectory("vendor/ddc/" "ddc") # SYSTEM)
endif()

set(koliop_BUILD_TESTING OFF CACHE BOOL "")
set(koliop_ENABLE_Kokkos        "SUBDIRECTORY" CACHE STRING "")
set(koliop_ENABLE_KokkosKernels "SUBDIRECTORY" CACHE STRING "")
add_subdirectory("vendor/koliop/" "koliop") # SYSTEM)

###############################################################################################
#                              Build libraries and executables
###############################################################################################
# Our project code

## Change the C++ flags to the Gyselalib++ specific ones
set(CMAKE_CXX_FLAGS "${GYSELALIBXX_DEFAULT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set(GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx xperiod_vx_user_defined_grid xnonperiod_vx_user_defined_grid)

## Get the main code in `src/`
add_subdirectory(src/)

## Build the simulations (applications mains) in `simulations/`
if("${GYSELALIBXX_BUILD_SIMULATIONS}")
add_subdirectory(simulations/)
endif()

## if tests are enabled, build the tests in `tests/`
if("${GYSELALIBXX_BUILD_TESTING}")
    add_subdirectory(tests/)
endif()

endif() # GYSELALIBXX_COMPILE_SOURCE
