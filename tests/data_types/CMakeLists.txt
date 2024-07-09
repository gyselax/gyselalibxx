# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_data_types
    deriv_field.cpp
    device_host_t.cpp
    field.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_data_types
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::data_types
)

gtest_discover_tests(unit_tests_data_types)
