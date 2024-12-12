# SPDX-License-Identifier: MIT

include(GoogleTest)

# Test advection for the operator in bsl_advection_1d.hpp
add_executable(1d_advection
    1d_advection_x.cpp
    1d_advection_xvx.cpp
    1d_advection_xyvxvy.cpp
    spatial_advection_1d.cpp
    velocity_advection_1d.cpp
    ../main.cpp
)
target_link_libraries(1d_advection
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        sll::SLL

        gslx::advection
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
)
target_compile_definitions(1d_advection PUBLIC)


gtest_discover_tests(1d_advection
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)
