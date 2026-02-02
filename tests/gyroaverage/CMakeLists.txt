# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_gyroaverage
    gyroaverage_circular.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_gyroaverage
    PUBLIC
        GTest::gtest
        GTest::gmock
        DDC::core
        DDC::splines
        gslx::utils
        gslx::gyroaverage
        testing_tools
)
gtest_discover_tests(unit_tests_gyroaverage DISCOVERY_MODE PRE_TEST)

