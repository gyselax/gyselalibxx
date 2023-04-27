# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.15)

add_executable(ddc_tests
    main.cpp
    aligned_allocator.cpp
    chunk.cpp
    discrete_domain.cpp
    non_uniform_point_sampling.cpp
    single_discretization.cpp
    tagged_vector.cpp
    relocatable_device_code_initialization.cpp
    relocatable_device_code.cpp
    type_seq.cpp
    uniform_point_sampling.cpp
    transform_reduce.cpp
    rectilinear_domain.cpp
    uniform_domain.cpp
    for_each.cpp
    fill.cpp
    discrete_element.cpp
    discrete_vector.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
    PUBLIC
        GTest::gtest
        DDC::DDC
)
gtest_discover_tests(ddc_tests)
