# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_io
    test_pdi.cpp
)
target_link_libraries(unit_tests_io
    PUBLIC
        DDC::pdi
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        PDI::pdi
        gslx::data_types
        gslx::io
        gslx::utils
)

gtest_discover_tests(unit_tests_io DISCOVERY_MODE PRE_TEST)
