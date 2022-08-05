# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)
project(voicexx C CXX)

# Our CMake module path

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")

# Configuration

## We use C++17 all around!
set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "The C++ standard whose features are requested to build this project.")

## List of options
option(BUILD_BENCHMARKS "Build the benchmarks." OFF)
option(VOICEXX_ENABLE_DEPRECATED "Enable deprecated code" OFF)
set(VOICEXX_DEFAULT_CXX_FLAGS "-O1" CACHE STRING "Default flags for C++ specific to Voice++")

# Our dependencies

## Use CTest for running tests
include(CTest)

## Look for a pre-installed LAPACK
find_package(LAPACK REQUIRED COMPONENTS CXX)

## Look for a pre-installed FFTW
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)

## Look for a pre-installed PDI
find_package(PDI REQUIRED COMPONENTS C)

## if tests are enabled, use googletest from `vendor/`
### we use it to write unit tests
if("${BUILD_TESTING}")
    option(INSTALL_GTEST "Enable installation of googletest." OFF)
    add_subdirectory("vendor/googletest/" "googletest")
endif()

## if benchmarks are enabled, use (google) benchmark from `vendor/`
### we use it to write our micro-benchmarks
### this must come after googletest
if(${BUILD_BENCHMARKS})
    option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
    add_subdirectory("vendor/benchmark/" "benchmark")
endif()

## Use mdspan from `vendor/`
set(MDSPAN_CXX_STANDARD "${CMAKE_CXX_STANDARD}" CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
add_subdirectory("vendor/mdspan/" "mdspan")

## Use Kokkos from `vendor/`
add_subdirectory("vendor/kokkos/" "kokkos")

## Use the discrete domain computation library (ddc) from `vendor/`
add_subdirectory("vendor/ddc/" "ddc")

## Use selalib's splines from `vendor/`
add_subdirectory("vendor/sll/" "sll")


# Our project code

## Change the C++ flags to the voice specific ones
set(CMAKE_CXX_FLAGS "${VOICEXX_DEFAULT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

## Get the main code in `src/`
set(GEOMETRY_VARIANTS_LIST periodic_rdimx nonperiodic_rdimx)
add_subdirectory(src/)

## Build the simulations (applications mains) in `simulations/`
add_subdirectory(simulations/landau)
add_subdirectory(simulations/bump_on_tail)

## if tests are enabled, build the tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
