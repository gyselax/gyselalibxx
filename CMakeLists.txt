cmake_minimum_required(VERSION 3.15)
project(voicexx CXX)


# Configuration

## We use C++17 all around, force it!
set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "" FORCE)


# Our dependencies

## Use CTest for running tests
include(CTest)

## Look for a pre-installed lapack
find_package(LAPACK REQUIRED COMPONENTS CXX)

## Look for googletest in `vendor/`, we use it to write unit tests
option(INSTALL_GTEST "Enable installation of googletest." OFF)
add_subdirectory("vendor/googletest/" "googletest")

## Look for mdspan in `vendor/`
set(MDSPAN_CXX_STANDARD 17 CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
add_subdirectory("vendor/mdspan/" "mdspan")

## Look for mdarray in `vendor/`
set(MDARRAY_CXX_STANDARD 17 CACHE STRING "Override the default CXX_STANDARD to compile with." FORCE)
add_subdirectory("vendor/mdarray/" "mdarray")


# Our project code

## The main code in `src/`
add_subdirectory(src)

## if tests are enabled, tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
