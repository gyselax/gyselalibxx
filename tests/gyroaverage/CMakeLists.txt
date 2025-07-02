# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_gyroaverage
    gyroaverage_circular.cpp
    ../main.cpp
)
target_compile_features(unit_tests_gyroaverage PUBLIC cxx_std_17)
target_link_libraries(unit_tests_gyroaverage
    PUBLIC
        GTest::gtest
        GTest::gmock
        DDC::core
        DDC::splines
        gslx::utils
        gslx::gyroaverage
)
gtest_discover_tests(unit_tests_gyroaverage DISCOVERY_MODE PRE_TEST)

