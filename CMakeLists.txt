# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.20)
project(gyselalibxx C CXX)

set_property(GLOBAL PROPERTY GLOBAL_DEPENDS_NO_CYCLES ON)

# Our CMake module path

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")

# Configuration

## We use C++17 all around!
set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "The C++ standard whose features are requested to build this project.")

## List of options
option(BUILD_BENCHMARKS "Build the benchmarks." OFF)
option(BUILD_DOCUMENTATION "Build the documentation." OFF)
option(GYSELALIBXX_ENABLE_DEPRECATED "Enable deprecated code" OFF)
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "-O1" CACHE STRING "Default flags for C++ specific to Voice++")

set(GYSELALIBXX_DEPENDENCY_POLICIES "AUTO" "EMBEDDED" "INSTALLED")

# Our dependencies

## Use CTest for running tests
include(CTest)

## Look for a pre-installed LAPACK
find_package(LAPACK REQUIRED COMPONENTS CXX)

## Look for a pre-installed paraconf
find_package(paraconf REQUIRED COMPONENTS C)

## Look for a pre-installed PDI
find_package(PDI REQUIRED COMPONENTS C)

## if tests are enabled, use googletest from `vendor/`
### we use it to write unit tests
if("${BUILD_TESTING}")
  set(GYSELALIBXX_GTest_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `GTest` package. Options: ${GYSELALIBXX_DEPENDENCY_POLICIES}")
  set_property(CACHE GYSELALIBXX_GTest_DEPENDENCY_POLICY PROPERTY STRINGS ${GYSELALIBXX_DEPENDENCY_POLICIES})

  if("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "AUTO")
    find_package(GTest "1.12" QUIET)
    if(NOT "${GTest_FOUND}")
      add_subdirectory("vendor/googletest/" "googletest")
    endif()
  elseif("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
    add_subdirectory("vendor/googletest/" "googletest")
  elseif("${GYSELALIBXX_GTest_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
    find_package(GTest "1.12" REQUIRED)
  endif()
endif()

## if benchmarks are enabled, use (google) benchmark from `vendor/`
### we use it to write our micro-benchmarks
### this must come after googletest
if("${BUILD_BENCHMARKS}")
    set(GYSELALIBXX_benchmark_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `benchmark` package. Options: ${GYSELALIBXX_DEPENDENCY_POLICIES}")
    set_property(CACHE GYSELALIBXX_benchmark_DEPENDENCY_POLICY PROPERTY STRINGS ${GYSELALIBXX_DEPENDENCY_POLICIES})

    if("${GYSELALIBXX_benchmark_DEPENDENCY_POLICY}" STREQUAL "AUTO")
      find_package(benchmark "1.6" QUIET)
      if(NOT "${benchmark_FOUND}")
        option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
        option(BENCHMARK_ENABLE_INSTALL "Enable installation of benchmark. (Projects embedding benchmark may want to turn this OFF.)" OFF)
        add_subdirectory("vendor/benchmark/" "benchmark")
      endif()
    elseif("${GYSELALIBXX_benchmark_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
      option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
      option(BENCHMARK_ENABLE_INSTALL "Enable installation of benchmark. (Projects embedding benchmark may want to turn this OFF.)" OFF)
      add_subdirectory("vendor/benchmark/" "benchmark")
    elseif("${GYSELALIBXX_benchmark_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
      find_package(benchmark "1.6" REQUIRED)
    endif()
endif()

## Use mdspan from `vendor/`
set(GYSELALIBXX_mdspan_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `mdspan` package. Options: ${GYSELALIBXX_DEPENDENCY_POLICIES}")
set_property(CACHE GYSELALIBXX_mdspan_DEPENDENCY_POLICY PROPERTY STRINGS ${GYSELALIBXX_DEPENDENCY_POLICIES})

if("${GYSELALIBXX_mdspan_DEPENDENCY_POLICY}" STREQUAL "AUTO")
  find_package(mdspan "0.3" QUIET)
  if(NOT "${mdspan_FOUND}")
    set(MDSPAN_CXX_STANDARD "${CMAKE_CXX_STANDARD}" CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
    add_subdirectory("vendor/mdspan/" "mdspan")
  endif()
elseif("${GYSELALIBXX_mdspan_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
  set(MDSPAN_CXX_STANDARD "${CMAKE_CXX_STANDARD}" CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
  add_subdirectory("vendor/mdspan/" "mdspan")
elseif("${GYSELALIBXX_mdspan_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
  find_package(mdspan "0.3" REQUIRED)
endif()

## if documentation is enabled, build it
if("${BUILD_DOCUMENTATION}")
  ## Look for a pre-installed Doxygen
  find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

  add_subdirectory(docs/)
endif()

## Turn off documentation for subpackages
set(BUILD_DOCUMENTATION OFF)

## Use Kokkos from `vendor/`
add_subdirectory("vendor/kokkos/" "kokkos")

find_package (Eigen3 3.3 NO_MODULE)

if (NOT EIGEN3_FOUND)
  ## Use Eigen from `vendor/`
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "Don't build eigen docs")
  add_subdirectory("vendor/eigen/" "eigen")
endif (NOT EIGEN3_FOUND)

## Use the discrete domain computation library (ddc) from `vendor/`
add_subdirectory("vendor/ddc/" "ddc")

## Use selalib's splines from `vendor/`
add_subdirectory("vendor/sll/" "sll")


# Our project code

## Change the C++ flags to the voice specific ones
set(CMAKE_CXX_FLAGS "${GYSELALIBXX_DEFAULT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set(GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

## Get the main code in `src/`
add_subdirectory(src/)

## Build the simulations (applications mains) in `simulations/`
add_subdirectory(simulations/)

## if tests are enabled, build the tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
