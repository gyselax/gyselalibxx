# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.15)

add_executable(ddc_tests
    aligned_allocator.cpp
    chunk.cpp
    discrete_domain.cpp
    non_uniform_discretization.cpp
    single_discretization.cpp
    tagged_vector.cpp
    type_seq.cpp
    uniform_discretization.cpp
    transform_reduce.cpp
    rectilinear_domain.cpp
    uniform_domain.cpp
    for_each.cpp
    discrete_element.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
    PUBLIC
        GTest::gtest_main
        DDC::DDC
)
gtest_discover_tests(ddc_tests)
