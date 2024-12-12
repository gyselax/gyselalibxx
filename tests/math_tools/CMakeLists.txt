# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_math_tools
    test_lnorm_tools.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_math_tools
    PUBLIC
        DDC::pdi
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
	gslx::math_tools
        gslx::utils
)

gtest_discover_tests(unit_tests_math_tools DISCOVERY_MODE PRE_TEST)
