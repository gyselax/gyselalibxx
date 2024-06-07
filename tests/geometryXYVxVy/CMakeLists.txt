# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

add_executable(unit_tests_xy_vxvy
    ../main.cpp
    quadrature.cpp
)
target_compile_features(unit_tests_xy_vxvy PUBLIC cxx_std_17)
target_link_libraries(unit_tests_xy_vxvy
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        gslx::geometry_xyvxvy
        sll::SLL
        gslx::quadrature
)

gtest_discover_tests(unit_tests_xy_vxvy
    TEST_SUFFIX "_xy_vxvy"
    DISCOVERY_MODE PRE_TEST)


add_subdirectory(landau)
