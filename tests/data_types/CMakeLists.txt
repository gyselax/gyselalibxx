# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests_data_types
    deriv_field.cpp
    device_host_t.cpp
    field.cpp
    ../main.cpp
)
target_compile_features(unit_tests_data_types PUBLIC cxx_std_17)
target_link_libraries(unit_tests_data_types
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::data_types
)

gtest_discover_tests(unit_tests_data_types)
