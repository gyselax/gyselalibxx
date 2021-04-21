cmake_minimum_required(VERSION 3.15)

project(voicexx CXX)

set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "" FORCE)

include(CTest)

find_package(LAPACK REQUIRED COMPONENTS CXX)

set(MDSPAN_CXX_STANDARD 17 CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
add_subdirectory("vendor/mdspan/" "mdspan" EXCLUDE_FROM_ALL)
set(MDARRAY_CXX_STANDARD 17 CACHE STRING "Override the default CXX_STANDARD to compile with" FORCE)
add_subdirectory("vendor/mdarray/" "mdarray" EXCLUDE_FROM_ALL)

add_subdirectory(src)
if("${BUILD_TESTING}")
add_subdirectory(tests/)
endif()
