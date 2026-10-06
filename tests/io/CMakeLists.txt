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
        PDI::PDI_C
        gyselalibxx::data_types
        gyselalibxx::io
        gyselalibxx::utils
)

gtest_discover_tests(unit_tests_io DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
