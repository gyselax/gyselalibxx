# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_utils
    test_ddcHelpers.cpp
    test_math_tools.cpp
    transpose.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_utils
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
	gslx::math_tools
        gslx::utils
)

gtest_discover_tests(unit_tests_utils DISCOVERY_MODE PRE_TEST)
