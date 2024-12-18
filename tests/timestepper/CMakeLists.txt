# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_timestepper
    euler_1d.cpp
    crank_nicolson_1d.cpp
    runge_kutta_1d.cpp
    runge_kutta_2d.cpp
    runge_kutta_2d_mixed.cpp
    euler_2d_mixed.cpp
    crank_nicolson_2d_mixed.cpp
    timestepper_multipatch.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_timestepper
    PUBLIC
        DDC::pdi
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        gslx::multipatch_data_types
        gslx::multipatch_geometries
        gslx::timestepper
        gslx::utils
)

gtest_discover_tests(unit_tests_timestepper DISCOVERY_MODE PRE_TEST)
