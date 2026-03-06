# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_utils
    non_uniform_interpolation_points.cpp
    spline_builder_deriv_field_2d_test.cpp
    test_ddcHelpers.cpp
    transpose.cpp
    view.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_utils
    PUBLIC
        DDC::pdi
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        gslx::utils
        gslx::data_types
)

gtest_discover_tests(unit_tests_utils DISCOVERY_MODE PRE_TEST)
target_include_directories(unit_tests_utils
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)
