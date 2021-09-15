cmake_minimum_required(VERSION 3.15)
project(voicexx CXX)

# Our include path

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")


# Configuration

## We use C++17 all around!
set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "The C++ standard whose features are requested to build this project.")

# List of options
option(BUILD_BENCHMARKS "Build the benchmarks." OFF)

# Our dependencies

## Use CTest for running tests
include(CTest)

## Look for a pre-installed lapack
find_package(LAPACK REQUIRED COMPONENTS CXX)

## Look for FFTW
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)

find_package(paraconf REQUIRED COMPONENTS C)

find_package(PDI REQUIRED COMPONENTS C)

## if tests are enabled, look for googletest in `vendor/`, we use it to write unit tests
if("${BUILD_TESTING}")
    option(INSTALL_GTEST "Enable installation of googletest." OFF)
    add_subdirectory("vendor/googletest/" "googletest")
endif()

## Must come after googletest
if(${BUILD_BENCHMARKS})
    add_subdirectory("vendor/benchmark/" "benchmark")
endif()

## Look for mdspan in `vendor/`
set(MDSPAN_CXX_STANDARD "${CMAKE_CXX_STANDARD}" CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
add_subdirectory("vendor/mdspan/" "mdspan")

## Look for the Domain DeComposition library (ddc) in `vendor/`
add_subdirectory("vendor/ddc/" "ddc")

## Look for selalib's splines in `vendor/`
add_subdirectory("vendor/sll/" "sll")


# Our project code

## The main code in `src/`
add_subdirectory(src)

## if tests are enabled, build the tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
