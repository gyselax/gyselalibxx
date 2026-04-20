# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)

# Read version from external file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/src/version.txt" GYSELALIBXX_VERSION_STRING)
string(STRIP "${GYSELALIBXX_VERSION_STRING}" GYSELALIBXX_VERSION_STRING)

# Set project language and version
project(gyselalibxx
    LANGUAGES C CXX
    VERSION ${GYSELALIBXX_VERSION_STRING})

set_property(GLOBAL PROPERTY GLOBAL_DEPENDS_NO_CYCLES ON)

# Configuration

# Ensure custom cmake modules can be found
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

## We use C++20 all around!
set(CMAKE_CXX_STANDARD 20 CACHE INTERNAL "The C++ standard whose features are requested to build this project.")

## List of options
option(GYSELALIBXX_BUILD_SIMULATIONS "Build the simulations." ON)
option(GYSELALIBXX_BUILD_TESTING "Build the tests." ON)
option(GYSELALIBXX_ENABLE_DEPRECATED "Enable deprecated code" OFF)
option(GYSELALIBXX_COMPILE_SOURCE "Enable compilation of the source code (this should be set to off to build documentation without the C++ dependencies)" ON)
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "-O1" CACHE STRING "Default flags for C++ specific to Gyselalib++")
option(GYSELALIBXX_ACTIVATE_RESTART_TESTS "Activate tests which check that a simulation gives the same results after restart." ON)

###############################################################################################
#                                     Get dependencies
###############################################################################################

if ("${GYSELALIBXX_COMPILE_SOURCE}")

# Our dependencies

## Use CTest for running tests
include(CTest)

find_package(LAPACK REQUIRED)

find_package(MPI REQUIRED)

find_package(paraconf 1.0.0...<2 REQUIRED COMPONENTS C)

find_package(PDI 1.10.1...<2 REQUIRED COMPONENTS C)

find_package(GTest 1.12 REQUIRED)

find_package(Kokkos 4.4.1...<5 REQUIRED)

find_package(KokkosKernels 4.5.1...<5 REQUIRED)

find_package(DDC 0.12.0 REQUIRED COMPONENTS fft pdi splines)

find_package(koliop 0.1.2 REQUIRED)

find_package(LAPACKE REQUIRED)

find_package(Ginkgo 1.8 REQUIRED)

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
