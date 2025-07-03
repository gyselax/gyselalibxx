# SPDX-License-Identifier: MIT

include(GoogleTest)

# Test advection for the operator in bsl_advection_1d.hpp
add_executable(unit_tests_advection
    1d_advection_x.cpp
    1d_advection_xvx.cpp
    1d_advection_xyvxvy.cpp
    polar_foot_finder.cpp
    spatial_advection_1d.cpp
    velocity_advection_1d.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_advection
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gslx::advection
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
)
target_compile_definitions(unit_tests_advection PUBLIC)


gtest_discover_tests(unit_tests_advection
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)
