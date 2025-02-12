# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_xperiod_vx
    krooksource.cpp
    maxwellian.cpp
    ../main.cpp
)

target_link_libraries(unit_tests_xperiod_vx
    PUBLIC
        GTest::gtest
        gslx::quadrature
        gslx::rhs_xperiod_vx
        gslx::speciesinfo
        gslx::utils_xperiod_vx
)

gtest_discover_tests(unit_tests_xperiod_vx
    TEST_SUFFIX "_xperiod_vx"
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)
