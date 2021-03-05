cmake_minimum_required(VERSION 3.15)

project(voicexx CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
add_subdirectory("vendor/mdspan/" "mdspan" EXCLUDE_FROM_ALL)

add_subdirectory(src)
