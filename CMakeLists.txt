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
option(GYSELALIBXX_BUILD_BENCHMARKS "Build the benchmarks." OFF)
option(GYSELALIBXX_BUILD_DOCUMENTATION "Build the documentation." OFF)
option(GYSELALIBXX_BUILD_SIMULATIONS "Build the simulations." ON)
option(GYSELALIBXX_BUILD_TESTING "Build the tests." ON)
option(GYSELALIBXX_ENABLE_DEPRECATED "Enable deprecated code" OFF)
option(GYSELALIBXX_COMPILE_SOURCE "Enable compilation of the source code (this should be set to off to build documentation without the C++ dependencies)" ON)
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "-O1" CACHE STRING "Default flags for C++ specific to Voice++")
option(ACTIVATE_RESTART_TESTS "Activate tests which check that a simulation gives the same results after restart." ON)

set(GYSELALIBXX_DEPENDENCY_POLICIES "AUTO" "EMBEDDED" "INSTALLED")

# Set default DDC options when included in GYSELA
option(DDC_BUILD_BENCHMARKS       "Build DDC benchmarks." OFF)
option(DDC_BUILD_DOCUMENTATION    "Build DDC documentation/website" OFF)
option(DDC_BUILD_EXAMPLES         "Build DDC examples" OFF)
option(DDC_BUILD_TESTS            "Build DDC tests if BUILD_TESTING is enabled" OFF)

###############################################################################################
#                                  Build documentation
###############################################################################################

## if documentation is enabled, build it
if("${GYSELALIBXX_BUILD_DOCUMENTATION}")
  ## Look for a pre-installed Doxygen
  find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

  add_subdirectory(docs/)
endif()

## Turn off documentation for subpackages
set(BUILD_DOCUMENTATION OFF)

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
if("${GYSELALIBXX_BUILD_TESTING}")
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

## if benchmarks are enabled, use (google) benchmark from `vendor/`
### we use it to write our micro-benchmarks
### this must come after googletest
if("${GYSELALIBXX_BUILD_BENCHMARKS}")
    set(GYSELALIBXX_benchmark_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `benchmark` package. Options: ${GYSELALIBXX_DEPENDENCY_POLICIES}")
    set_property(CACHE GYSELALIBXX_benchmark_DEPENDENCY_POLICY PROPERTY STRINGS ${GYSELALIBXX_DEPENDENCY_POLICIES})

    if("${GYSELALIBXX_benchmark_DEPENDENCY_POLICY}" STREQUAL "AUTO")
      find_package(benchmark "1.6" QUIET)
      if(NOT "${benchmark_FOUND}")
        option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
        option(BENCHMARK_ENABLE_INSTALL "Enable installation of benchmark. (Projects embedding benchmark may want to turn this OFF.)" OFF)
        add_subdirectory("vendor/benchmark/" "benchmark") # SYSTEM)
      endif()
    elseif("${GYSELALIBXX_benchmark_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
      option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
      option(BENCHMARK_ENABLE_INSTALL "Enable installation of benchmark. (Projects embedding benchmark may want to turn this OFF.)" OFF)
      add_subdirectory("vendor/benchmark/" "benchmark") # SYSTEM)
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
    add_subdirectory("vendor/mdspan/" "mdspan") # SYSTEM)
  endif()
elseif("${GYSELALIBXX_mdspan_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
  set(MDSPAN_CXX_STANDARD "${CMAKE_CXX_STANDARD}" CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
  add_subdirectory("vendor/mdspan/" "mdspan") # SYSTEM)
elseif("${GYSELALIBXX_mdspan_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
  find_package(mdspan "0.3" REQUIRED)
endif()

## Use Kokkos from `vendor/`
if("${Kokkos_ENABLE_CUDA}")
  option(Kokkos_ENABLE_CUDA_CONSTEXPR "Whether to activate experimental relaxed constexpr functions" ON)
  option(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for CUDA" ON)
endif()
add_subdirectory("vendor/kokkos/" "kokkos") # SYSTEM)

# NOTE: Do not forget to keep Kokkos and Kokkos-kernels in sync !
set(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "")
add_subdirectory("vendor/kokkos-kernels/" "kokkos-kernels") # SYSTEM)

find_package (Eigen3 3.3 NO_MODULE)

if (NOT EIGEN3_FOUND)
  ## Use Eigen from `vendor/`
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "Don't build eigen docs")
  add_subdirectory("vendor/eigen/" "eigen") # SYSTEM)
endif (NOT EIGEN3_FOUND)

## Use the discrete domain computation library (ddc) from `vendor/`
add_subdirectory("vendor/ddc/" "ddc") # SYSTEM)

## Use selalib's splines from `vendor/`
add_subdirectory("vendor/sll/" "sll") # SYSTEM)

set(koliop_BUILD_TESTING OFF CACHE BOOL "")
set(koliop_ENABLE_Kokkos        "SUBDIRECTORY" CACHE STRING "")
set(koliop_ENABLE_KokkosKernels "SUBDIRECTORY" CACHE STRING "")
add_subdirectory("vendor/koliop/" "koliop") # SYSTEM)

## Use the kokkos-tools for profiling purpose from `vendor/`
include(vendor/kokkos-tools/cmake/utils.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vendor/kokkos-tools/common ${CMAKE_CURRENT_SOURCE_DIR}/vendor/kokkos-tools/profiling/all)
configure_file(vendor/kokkos-tools/common/kp_config.hpp.in kokkos-tools/profiling/simple-kernel-timer/kp_config.hpp)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/kokkos-tools/profiling/simple-kernel-timer) # Weird trick to get access to kp_config.hpp
set(BUILD_SHARED_LIBS_DUMP "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS ON) #Needed to get libkp_kernel_timer.so
add_subdirectory("vendor/kokkos-tools/profiling/simple-kernel-timer/" "kokkos-tools/profiling/simple-kernel-timer")
add_compile_definitions(KP_KERNEL_TIMER_PATH="${CMAKE_CURRENT_BINARY_DIR}/kokkos-tools/profiling/simple-kernel-timer/libkp_kernel_timer.so")
if (${KOKKOS_ENABLE_CUDA})
    add_subdirectory("vendor/kokkos-tools/profiling/nvtx-connector/" "kokkos-tools/profiling/nvtx-connector")
    add_compile_definitions(KOKKOS_TOOLS_LIBS="${CMAKE_CURRENT_BINARY_DIR}/kokkos-tools/profiling/nvtx-connector/libkp_nvtx_connector.so")
endif()
set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_DUMP}")

###############################################################################################
#                              Build libraries and executables
###############################################################################################
# Our project code

## Change the C++ flags to the voice specific ones
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
