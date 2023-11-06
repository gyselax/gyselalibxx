# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests_timestepper
	euler_1d.cpp
	crank_nicolson_1d.cpp
    runge_kutta_1d.cpp
    runge_kutta_2d.cpp
    runge_kutta_2d_mixed.cpp
    euler_2d_mixed.cpp
    crank_nicolson_2d_mixed.cpp
    ../main.cpp
)
target_compile_features(unit_tests_timestepper PUBLIC cxx_std_17)
target_link_libraries(unit_tests_timestepper
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        vcx::timestepper
        vcx::utils
)

gtest_discover_tests(unit_tests_timestepper)
