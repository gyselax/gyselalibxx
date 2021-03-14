cmake_minimum_required(VERSION 3.15)

project(voicexx CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(LAPACK REQUIRED COMPONENTS CXX)
add_subdirectory("vendor/mdspan/" "mdspan" EXCLUDE_FROM_ALL)

add_subdirectory(src)
