# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests_utils
    field.cpp
    ../main.cpp
)
target_compile_features(unit_tests_utils PUBLIC cxx_std_17)
target_link_libraries(unit_tests_utils
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        gslx::utils
)

gtest_discover_tests(unit_tests_utils)
