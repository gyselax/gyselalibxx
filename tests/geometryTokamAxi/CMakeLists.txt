# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_geometryTokamAxi
    maxwellian.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_geometryTokamAxi
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::geometry_tokamaxi
        gslx::initialization_tokamaxi
        gslx::quadrature
        gslx::utils
)

gtest_discover_tests(unit_tests_geometryTokamAxi DISCOVERY_MODE PRE_TEST)
